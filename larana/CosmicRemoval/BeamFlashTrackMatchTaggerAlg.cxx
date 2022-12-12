/*!
 * Title:   Beam Flash<-->Track Match Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code from Ben Jones
 *
 * Description: Algorithm that compares all tracks to the flash during the
 *              beam gate, and determines if that track is consistent with
 *              having produced that flash.
 * Input:       recob::OpFlash, recob::Track
 * Output:      anab::CosmicTag (and Assn<anab::CosmicTag,recob::Track>)
*/

#include "BeamFlashTrackMatchTaggerAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include "larana/OpticalDetector/OpDigiProperties.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include <limits>
#include <utility>
#include <vector>

#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"

cosmic::BeamFlashTrackMatchTaggerAlg::BeamFlashTrackMatchTaggerAlg(fhicl::ParameterSet const& p)
  : COSMIC_TYPE_FLASHMATCH(anab::CosmicTagID_t::kFlash_BeamIncompatible)
  , COSMIC_TYPE_OUTSIDEDRIFT(anab::CosmicTagID_t::kOutsideDrift_Partial)
  , DEBUG_FLAG(p.get<bool>("RunDebugMode", false))
  , fMinTrackLength(p.get<float>("MinTrackLength"))
  , fMinOpHitPE(p.get<float>("MinOpHitPE", 0.1))
  , fMIPdQdx(p.get<float>("MIPdQdx", 2.1))
  , fOpDetSaturation(p.get<float>("OpDetSaturation", 200.))
  , fSingleChannelCut(p.get<float>("SingleChannelCut"))
  , fCumulativeChannelThreshold(p.get<float>("CumulativeChannelThreshold"))
  , fCumulativeChannelCut(p.get<unsigned int>("CumulativeChannelCut"))
  , fIntegralCut(p.get<float>("IntegralCut"))
  , fMakeOutsideDriftTags(p.get<bool>("MakeOutsideDriftTags", false))
  , fNormalizeHypothesisToFlash(p.get<bool>("NormalizeHypothesisToFlash"))
{}

void cosmic::BeamFlashTrackMatchTaggerAlg::SetHypothesisComparisonTree(TTree* tree,
                                                                       TH1F* hist_flash,
                                                                       TH1F* hist_hyp)
{
  cTree = tree;

  cOpDetHist_flash = hist_flash;
  cOpDetHist_flash->SetNameTitle("opdet_hist_flash", "Optical Detector Occupancy, Flash");

  cOpDetHist_hyp = hist_hyp;
  cOpDetHist_hyp->SetNameTitle("opdet_hist_hyp", "Optical Detector Occupancy, Hypothesis");

  cTree->Branch("fcp", &cFlashComparison_p, cFlashComparison_p.leaf_structure.c_str());
  cTree->Branch("opdet_hyp", &cOpDetVector_hyp);
  cTree->Branch("opdet_flash", &cOpDetVector_flash);
  cTree->Branch("opdet_hist_flash", "TH1F", cOpDetHist_flash);
  cTree->Branch("opdet_hist_hyp", "TH1F", cOpDetHist_hyp);
}

void cosmic::BeamFlashTrackMatchTaggerAlg::RunCompatibilityCheck(
  std::vector<recob::OpFlash> const& flashVector,
  std::vector<recob::Track> const& trackVector,
  std::vector<anab::CosmicTag>& cosmicTagVector,
  std::vector<size_t>& assnTrackTagVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip)
{
  auto const& geom = *providers.get<geo::GeometryCore>();
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

  std::vector<const recob::OpFlash*> flashesOnBeamTime;
  for (auto const& flash : flashVector) {
    if (!flash.OnBeamTime()) continue;
    flashesOnBeamTime.push_back(&flash);
  }

  //make sure this association vector is initialized properly
  assnTrackTagVector.resize(trackVector.size(), std::numeric_limits<size_t>::max());
  cosmicTagVector.reserve(trackVector.size());

  for (size_t track_i = 0; track_i < trackVector.size(); track_i++) {

    recob::Track const& track(trackVector[track_i]);

    if (track.Length() < fMinTrackLength) continue;

    //get the begin and end points of this track
    TVector3 const& pt_begin = track.LocationAtPoint<TVector3>(0);
    TVector3 const& pt_end = track.LocationAtPoint<TVector3>(track.NumberTrajectoryPoints() - 1);
    std::vector<float> xyz_begin = {(float)pt_begin.x(), (float)pt_begin.y(), (float)pt_begin.z()};
    std::vector<float> xyz_end = {(float)pt_end.x(), (float)pt_end.y(), (float)pt_end.z()};

    //check if this track is outside the drift window, and if it is continue
    if (!InDriftWindow(pt_begin.x(), pt_end.x(), geom)) {
      if (fMakeOutsideDriftTags) {
        cosmicTagVector.emplace_back(xyz_begin, xyz_end, 1., COSMIC_TYPE_OUTSIDEDRIFT);
        assnTrackTagVector[track_i] = cosmicTagVector.size() - 1;
      }
      continue;
    }

    //get light hypothesis for track
    std::vector<float> lightHypothesis = GetMIPHypotheses(track, providers, pvs, opdigip);

    //check compatibility with beam flash
    bool compatible = false;
    for (const recob::OpFlash* flashPointer : flashesOnBeamTime) {
      CompatibilityResultType result =
        CheckCompatibility(lightHypothesis, flashPointer, geom, wireReadoutGeom);
      if (result == CompatibilityResultType::kCompatible) compatible = true;
      if (DEBUG_FLAG) {
        PrintTrackProperties(track);
        PrintFlashProperties(*flashPointer);
        PrintHypothesisFlashComparison(
          lightHypothesis, flashPointer, geom, wireReadoutGeom, result);
      }
    }

    //make tag
    float cosmicScore = 1.;
    if (compatible) cosmicScore = 0.;
    cosmicTagVector.emplace_back(xyz_begin, xyz_end, cosmicScore, COSMIC_TYPE_FLASHMATCH);
    assnTrackTagVector[track_i] = cosmicTagVector.size() - 1;
  }
}

//this compares the hypothesis to the flash itself.
void cosmic::BeamFlashTrackMatchTaggerAlg::RunHypothesisComparison(
  unsigned int const run,
  unsigned int const event,
  std::vector<recob::OpFlash> const& flashVector,
  std::vector<recob::Track> const& trackVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip)
{
  auto const& geom = *(providers.get<geo::GeometryCore>());
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

  cFlashComparison_p.run = run;
  cFlashComparison_p.event = event;

  std::vector<std::pair<unsigned int, const recob::OpFlash*>> flashesOnBeamTime;
  for (unsigned int i = 0; i < flashVector.size(); i++) {
    recob::OpFlash const& flash = flashVector[i];
    if (!flash.OnBeamTime()) continue;
    flashesOnBeamTime.push_back(std::make_pair(i, &flash));
  }

  for (size_t track_i = 0; track_i < trackVector.size(); track_i++) {

    recob::Track const& track(trackVector[track_i]);

    if (track.Length() < fMinTrackLength) continue;

    //get the begin and end points of this track
    TVector3 const& pt_begin = track.LocationAtPoint<TVector3>(0);
    TVector3 const& pt_end = track.LocationAtPoint<TVector3>(track.NumberTrajectoryPoints() - 1);
    std::vector<float> xyz_begin = {(float)pt_begin.x(), (float)pt_begin.y(), (float)pt_begin.z()};
    std::vector<float> xyz_end = {(float)pt_end.x(), (float)pt_end.y(), (float)pt_end.z()};

    //check if this track is outside the drift window, and if it is continue
    if (!InDriftWindow(pt_begin.x(), pt_end.x(), geom)) continue;

    cFlashComparison_p.trk_startx = pt_begin.x();
    cFlashComparison_p.trk_starty = pt_begin.y();
    cFlashComparison_p.trk_startz = pt_begin.z();
    cFlashComparison_p.trk_endx = pt_end.x();
    cFlashComparison_p.trk_endy = pt_end.y();
    cFlashComparison_p.trk_endz = pt_end.z();

    //get light hypothesis for track
    cOpDetVector_hyp = GetMIPHypotheses(track, providers, pvs, opdigip);

    cFlashComparison_p.hyp_index = track_i;
    FillFlashProperties(cOpDetVector_hyp,
                        cFlashComparison_p.hyp_totalPE,
                        cFlashComparison_p.hyp_y,
                        cFlashComparison_p.hyp_sigmay,
                        cFlashComparison_p.hyp_z,
                        cFlashComparison_p.hyp_sigmaz,
                        geom);

    for (auto flash : flashesOnBeamTime) {
      cOpDetVector_flash = std::vector<float>(geom.NOpDets(), 0);
      cFlashComparison_p.flash_nOpDet = 0;
      for (size_t c = 0; c <= wireReadoutGeom.MaxOpChannel(); c++) {
        if (wireReadoutGeom.IsValidOpChannel(c)) {
          unsigned int OpDet = wireReadoutGeom.OpDetFromOpChannel(c);
          cOpDetVector_flash[OpDet] += flash.second->PE(c);
        }
      }
      for (size_t o = 0; o < cOpDetVector_flash.size(); o++)
        if (cOpDetVector_flash[o] < fMinOpHitPE) cFlashComparison_p.flash_nOpDet++;

      cFlashComparison_p.flash_index = flash.first;
      cFlashComparison_p.flash_totalPE = flash.second->TotalPE();
      cFlashComparison_p.flash_y = flash.second->YCenter();
      cFlashComparison_p.flash_sigmay = flash.second->YWidth();
      cFlashComparison_p.flash_z = flash.second->ZCenter();
      cFlashComparison_p.flash_sigmaz = flash.second->ZWidth();

      cFlashComparison_p.chi2 = CalculateChi2(cOpDetVector_flash, cOpDetVector_hyp);

      cTree->Fill();
    } //end loop over flashes

  } //end loop over tracks
}

//this compares the hypothesis to the flash itself.
void cosmic::BeamFlashTrackMatchTaggerAlg::RunHypothesisComparison(
  unsigned int const run,
  unsigned int const event,
  std::vector<recob::OpFlash> const& flashVector,
  std::vector<simb::MCParticle> const& mcParticleVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip)
{
  auto const& geom = *(providers.get<geo::GeometryCore>());
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

  cFlashComparison_p.run = run;
  cFlashComparison_p.event = event;

  std::vector<std::pair<unsigned int, const recob::OpFlash*>> flashesOnBeamTime;
  for (unsigned int i = 0; i < flashVector.size(); i++) {
    recob::OpFlash const& flash = flashVector[i];
    if (!flash.OnBeamTime()) continue;
    flashesOnBeamTime.push_back(std::make_pair(i, &flash));
  }

  for (size_t particle_i = 0; particle_i < mcParticleVector.size(); particle_i++) {

    simb::MCParticle const& particle(mcParticleVector[particle_i]);
    if (particle.Process().compare("primary") != 0) continue;
    if (particle.Trajectory().TotalLength() < fMinTrackLength) continue;

    //get the begin and end points of this track
    size_t start_i = 0, end_i = particle.NumberTrajectoryPoints() - 1;
    bool prev_inside = false;
    for (size_t pt_i = 0; pt_i < particle.NumberTrajectoryPoints(); pt_i++) {
      bool inside = InDetector(particle.Position(pt_i).Vect(), geom);
      if (inside && !prev_inside) start_i = pt_i;
      if (!inside && prev_inside) {
        end_i = pt_i - 1;
        break;
      }
      prev_inside = inside;
    }
    TVector3 const& pt_begin = particle.Position(start_i).Vect();
    TVector3 const& pt_end = particle.Position(end_i).Vect();
    std::vector<float> xyz_begin = {(float)pt_begin.x(), (float)pt_begin.y(), (float)pt_begin.z()};
    std::vector<float> xyz_end = {(float)pt_end.x(), (float)pt_end.y(), (float)pt_end.z()};

    //check if this track is outside the drift window, and if it is continue
    if (!InDriftWindow(pt_begin.x(), pt_end.x(), geom)) continue;

    cFlashComparison_p.trk_startx = pt_begin.x();
    cFlashComparison_p.trk_starty = pt_begin.y();
    cFlashComparison_p.trk_startz = pt_begin.z();
    cFlashComparison_p.trk_endx = pt_end.x();
    cFlashComparison_p.trk_endy = pt_end.y();
    cFlashComparison_p.trk_endz = pt_end.z();

    //get light hypothesis for track
    cOpDetVector_hyp = GetMIPHypotheses(particle, start_i, end_i, providers, pvs, opdigip);

    cFlashComparison_p.hyp_index = particle_i;
    FillFlashProperties(cOpDetVector_hyp,
                        cFlashComparison_p.hyp_totalPE,
                        cFlashComparison_p.hyp_y,
                        cFlashComparison_p.hyp_sigmay,
                        cFlashComparison_p.hyp_z,
                        cFlashComparison_p.hyp_sigmaz,
                        geom);

    for (auto flash : flashesOnBeamTime) {
      cOpDetVector_flash = std::vector<float>(geom.NOpDets(), 0);
      cFlashComparison_p.flash_nOpDet = 0;
      for (size_t c = 0; c <= wireReadoutGeom.MaxOpChannel(); c++) {
        if (wireReadoutGeom.IsValidOpChannel(c)) {
          unsigned int OpDet = wireReadoutGeom.OpDetFromOpChannel(c);
          cOpDetVector_flash[OpDet] += flash.second->PE(c);
        }
      }
      for (size_t o = 0; o < cOpDetVector_flash.size(); o++)
        if (cOpDetVector_flash[o] < fMinOpHitPE) cFlashComparison_p.flash_nOpDet++;
      cFlashComparison_p.flash_index = flash.first;
      cFlashComparison_p.flash_totalPE = flash.second->TotalPE();
      cFlashComparison_p.flash_y = flash.second->YCenter();
      cFlashComparison_p.flash_sigmay = flash.second->YWidth();
      cFlashComparison_p.flash_z = flash.second->ZCenter();
      cFlashComparison_p.flash_sigmaz = flash.second->ZWidth();

      cFlashComparison_p.chi2 = CalculateChi2(cOpDetVector_flash, cOpDetVector_hyp);

      for (size_t i = 0; i < cOpDetVector_flash.size(); i++)
        cOpDetHist_flash->SetBinContent(i + 1, cOpDetVector_flash[i]);

      for (size_t i = 0; i < cOpDetVector_hyp.size(); i++)
        cOpDetHist_hyp->SetBinContent(i + 1, cOpDetVector_hyp[i]);

      for (size_t i = 0; i < cOpDetVector_flash.size(); i++) {
        std::cout << "Flash/Hyp " << i << " : " << cOpDetHist_flash->GetBinContent(i + 1) << " "
                  << cOpDetHist_hyp->GetBinContent(i + 1) << std::endl;
      }

      cTree->Fill();
    } //end loop over flashes

  } //end loop over tracks
}

void cosmic::BeamFlashTrackMatchTaggerAlg::FillFlashProperties(
  std::vector<float> const& opdetVector,
  float& sum,
  float& y,
  float& sigmay,
  float& z,
  float& sigmaz,
  geo::GeometryCore const& geom)
{
  y = 0;
  sigmay = 0;
  z = 0;
  sigmaz = 0;
  sum = 0;
  for (unsigned int opdet = 0; opdet < opdetVector.size(); opdet++) {
    sum += opdetVector[opdet];
    auto const xyz = geom.Cryostat().OpDet(opdet).GetCenter();
    y += opdetVector[opdet] * xyz.Y();
    z += opdetVector[opdet] * xyz.Z();
  }

  y /= sum;
  z /= sum;

  for (unsigned int opdet = 0; opdet < opdetVector.size(); opdet++) {
    auto const xyz = geom.Cryostat().OpDet(opdet).GetCenter();
    sigmay += (opdetVector[opdet] * xyz.Y() - y) * (opdetVector[opdet] * xyz.Y() - y);
    sigmaz += (opdetVector[opdet] * xyz.Z() - y) * (opdetVector[opdet] * xyz.Z() - y);
  }

  sigmay = std::sqrt(sigmay) / sum;
  sigmaz = std::sqrt(sigmaz) / sum;
}

bool cosmic::BeamFlashTrackMatchTaggerAlg::InDetector(TVector3 const& pt,
                                                      geo::GeometryCore const& geom)
{
  auto const& tpc = geom.TPC({0, 0});
  if (pt.x() < 0 || pt.x() > 2 * tpc.HalfWidth()) return false;
  if (std::abs(pt.y()) > tpc.HalfHeight()) return false;
  if (pt.z() < 0 || pt.z() > tpc.Length()) return false;
  return true;
}

bool cosmic::BeamFlashTrackMatchTaggerAlg::InDriftWindow(double start_x,
                                                         double end_x,
                                                         geo::GeometryCore const& geom)
{
  if (start_x < 0. || end_x < 0.) return false;
  auto const& tpc = geom.TPC({0, 0});
  if (start_x > 2 * tpc.HalfWidth() || end_x > 2 * tpc.HalfWidth()) return false;
  return true;
}

void cosmic::BeamFlashTrackMatchTaggerAlg::AddLightFromSegment(
  TVector3 const& pt1,
  TVector3 const& pt2,
  std::vector<float>& lightHypothesis,
  float& totalHypothesisPE,
  geo::GeometryCore const& geom,
  phot::PhotonVisibilityService const& pvs,
  float const& PromptMIPScintYield,
  float XOffset)
{
  double xyz_segment[3];
  xyz_segment[0] = 0.5 * (pt2.x() + pt1.x()) + XOffset;
  xyz_segment[1] = 0.5 * (pt2.y() + pt1.y());
  xyz_segment[2] = 0.5 * (pt2.z() + pt1.z());

  //get the visibility vector
  auto const& PointVisibility = pvs.GetAllVisibilities(xyz_segment);

  //check pointer, as it may be null if given a y/z outside some range
  if (!PointVisibility) return;

  //get the amount of light
  float LightAmount = PromptMIPScintYield * (pt2 - pt1).Mag();

  for (size_t opdet_i = 0; opdet_i < pvs.NOpChannels(); opdet_i++) {
    lightHypothesis[opdet_i] += PointVisibility[opdet_i] * LightAmount;
    totalHypothesisPE += PointVisibility[opdet_i] * LightAmount;

    //apply saturation limit
    if (lightHypothesis[opdet_i] > fOpDetSaturation) {
      totalHypothesisPE -= (lightHypothesis[opdet_i] - fOpDetSaturation);
      lightHypothesis[opdet_i] = fOpDetSaturation;
    }
  }

} //end AddLightFromSegment

void cosmic::BeamFlashTrackMatchTaggerAlg::NormalizeLightHypothesis(
  std::vector<float>& lightHypothesis,
  float const& totalHypothesisPE,
  geo::GeometryCore const& geom)
{
  for (size_t opdet_i = 0; opdet_i < geom.NOpDets(); opdet_i++)
    lightHypothesis[opdet_i] /= totalHypothesisPE;
}

// Get a hypothesis for the light collected for a track
std::vector<float> cosmic::BeamFlashTrackMatchTaggerAlg::GetMIPHypotheses(
  recob::Track const& track,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  auto const& geom = *(providers.get<geo::GeometryCore>());
  auto const& larp = *(providers.get<detinfo::LArProperties>());
  std::vector<float> lightHypothesis(geom.NOpDets(), 0);
  float totalHypothesisPE = 0;
  const float PromptMIPScintYield =
    larp.ScintYield() * larp.ScintYieldRatio() * opdigip.QE() * fMIPdQdx;

  //get QE from ubChannelConfig, which gives per tube, so goes in AddLightFromSegment
  //VisibleEnergySeparation(step);

  for (size_t pt = 1; pt < track.NumberTrajectoryPoints(); pt++)
    AddLightFromSegment(track.LocationAtPoint<TVector3>(pt - 1),
                        track.LocationAtPoint<TVector3>(pt),
                        lightHypothesis,
                        totalHypothesisPE,
                        geom,
                        pvs,
                        PromptMIPScintYield,
                        XOffset);

  if (fNormalizeHypothesisToFlash && totalHypothesisPE > std::numeric_limits<float>::epsilon())
    NormalizeLightHypothesis(lightHypothesis, totalHypothesisPE, geom);

  return lightHypothesis;

} //end GetMIPHypotheses

// Get a hypothesis for the light collected for a particle trajectory
std::vector<float> cosmic::BeamFlashTrackMatchTaggerAlg::GetMIPHypotheses(
  simb::MCParticle const& particle,
  size_t start_i,
  size_t end_i,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  auto const& geom = *(providers.get<geo::GeometryCore>());
  auto const& larp = *(providers.get<detinfo::LArProperties>());
  std::vector<float> lightHypothesis(geom.NOpDets(), 0);
  float totalHypothesisPE = 0;
  const float PromptMIPScintYield =
    larp.ScintYield() * larp.ScintYieldRatio() * opdigip.QE() * fMIPdQdx;

  for (size_t pt = start_i + 1; pt <= end_i; pt++)
    AddLightFromSegment(particle.Position(pt - 1).Vect(),
                        particle.Position(pt).Vect(),
                        lightHypothesis,
                        totalHypothesisPE,
                        geom,
                        pvs,
                        PromptMIPScintYield,
                        XOffset);

  if (fNormalizeHypothesisToFlash && totalHypothesisPE > std::numeric_limits<float>::epsilon())
    NormalizeLightHypothesis(lightHypothesis, totalHypothesisPE, geom);

  return lightHypothesis;

} //end GetMIPHypotheses

//---------------------------------------
//  Check whether a hypothesis can be accomodated in a flash
//   Flashes fail if 1 bin is far in excess of the observed signal
//   or if the whole flash intensity is much too large for the hypothesis.
//  MIP dEdx is assumed for now.  Accounting for real dQdx will
//   improve performance of this algorithm.
//---------------------------------------
cosmic::BeamFlashTrackMatchTaggerAlg::CompatibilityResultType
cosmic::BeamFlashTrackMatchTaggerAlg::CheckCompatibility(
  std::vector<float> const& lightHypothesis,
  const recob::OpFlash* flashPointer,
  geo::GeometryCore const& geom,
  geo::WireReadoutGeom const& wireReadoutGeom)
{
  float hypothesis_integral = 0;
  float flash_integral = 0;
  unsigned int cumulativeChannels = 0;

  std::vector<double> PEbyOpDet(geom.NOpDets(), 0);
  for (unsigned int c = 0; c <= wireReadoutGeom.MaxOpChannel(); c++) {
    if (wireReadoutGeom.IsValidOpChannel(c)) {
      unsigned int o = wireReadoutGeom.OpDetFromOpChannel(c);
      PEbyOpDet[o] += flashPointer->PE(c);
    }
  }

  float hypothesis_scale = 1.;
  if (fNormalizeHypothesisToFlash) hypothesis_scale = flashPointer->TotalPE();

  for (size_t pmt_i = 0; pmt_i < lightHypothesis.size(); pmt_i++) {

    flash_integral += PEbyOpDet[pmt_i];

    if (lightHypothesis[pmt_i] < std::numeric_limits<float>::epsilon()) continue;
    hypothesis_integral += lightHypothesis[pmt_i] * hypothesis_scale;

    if (PEbyOpDet[pmt_i] < fMinOpHitPE) continue;

    float diff_scaled = (lightHypothesis[pmt_i] * hypothesis_scale - PEbyOpDet[pmt_i]) /
                        std::sqrt(lightHypothesis[pmt_i] * hypothesis_scale);

    if (diff_scaled > fSingleChannelCut) return CompatibilityResultType::kSingleChannelCut;

    if (diff_scaled > fCumulativeChannelThreshold) cumulativeChannels++;
    if (cumulativeChannels >= fCumulativeChannelCut)
      return CompatibilityResultType::kCumulativeChannelCut;
  }

  if ((hypothesis_integral - flash_integral) / std::sqrt(hypothesis_integral) > fIntegralCut)
    return CompatibilityResultType::kIntegralCut;

  return CompatibilityResultType::kCompatible;
}

float cosmic::BeamFlashTrackMatchTaggerAlg::CalculateChi2(std::vector<float> const& light_flash,
                                                          std::vector<float> const& light_track)
{
  float chi2 = 0;
  for (size_t pmt_i = 0; pmt_i < light_flash.size(); pmt_i++) {

    if (light_flash[pmt_i] < fMinOpHitPE) continue;

    float err2 = 1;
    if (light_track[pmt_i] > 1) err2 = light_track[pmt_i];

    chi2 +=
      (light_flash[pmt_i] - light_track[pmt_i]) * (light_flash[pmt_i] - light_track[pmt_i]) / err2;
  }

  return chi2;
}

void cosmic::BeamFlashTrackMatchTaggerAlg::PrintTrackProperties(recob::Track const& track,
                                                                std::ostream* output)
{
  *output << "----------------------------------------------" << std::endl;
  *output << "Track properties: ";
  *output << "\n\tLength=" << track.Length();

  auto const& pt_begin = track.LocationAtPoint(0);
  *output << "\n\tBegin Location (x,y,z)=(" << pt_begin.x() << "," << pt_begin.y() << ","
          << pt_begin.z() << ")";

  auto const& pt_end = track.LocationAtPoint(track.NumberTrajectoryPoints() - 1);
  *output << "\n\tEnd Location (x,y,z)=(" << pt_end.x() << "," << pt_end.y() << "," << pt_end.z()
          << ")";

  *output << "\n\tTrajectoryPoints=" << track.NumberTrajectoryPoints();
  *output << std::endl;
  *output << "----------------------------------------------" << std::endl;
}

void cosmic::BeamFlashTrackMatchTaggerAlg::PrintFlashProperties(recob::OpFlash const& flash,
                                                                std::ostream* output)
{
  *output << "----------------------------------------------" << std::endl;
  *output << "Flash properties: ";

  *output << "\n\tTime=" << flash.Time();
  *output << "\n\tOnBeamTime=" << flash.OnBeamTime();
  *output << "\n\ty position (center,width)=(" << flash.YCenter() << "," << flash.YWidth() << ")";
  *output << "\n\tz position (center,width)=(" << flash.ZCenter() << "," << flash.ZWidth() << ")";
  *output << "\n\tTotal PE=" << flash.TotalPE();

  *output << std::endl;
  *output << "----------------------------------------------" << std::endl;
}

void cosmic::BeamFlashTrackMatchTaggerAlg::PrintHypothesisFlashComparison(
  std::vector<float> const& lightHypothesis,
  const recob::OpFlash* flashPointer,
  geo::GeometryCore const& geom,
  geo::WireReadoutGeom const& wireReadoutGeom,
  CompatibilityResultType result,
  std::ostream* output)
{
  *output << "----------------------------------------------" << std::endl;
  *output << "Hypothesis-flash comparison: ";

  float hypothesis_integral = 0;
  float flash_integral = 0;

  float hypothesis_scale = 1.;
  if (fNormalizeHypothesisToFlash) hypothesis_scale = flashPointer->TotalPE();

  std::vector<double> PEbyOpDet(geom.NOpDets(), 0);
  for (unsigned int c = 0; c <= wireReadoutGeom.MaxOpChannel(); c++) {
    if (wireReadoutGeom.IsValidOpChannel(c)) {
      unsigned int o = wireReadoutGeom.OpDetFromOpChannel(c);
      PEbyOpDet[o] += flashPointer->PE(c);
    }
  }

  for (size_t pmt_i = 0; pmt_i < lightHypothesis.size(); pmt_i++) {

    flash_integral += PEbyOpDet[pmt_i];

    *output << "\n\t pmt_i=" << pmt_i << ", (hypothesis,flash)=("
            << lightHypothesis[pmt_i] * hypothesis_scale << "," << PEbyOpDet[pmt_i] << ")";

    if (lightHypothesis[pmt_i] < std::numeric_limits<float>::epsilon()) continue;

    *output << "  difference="
            << (lightHypothesis[pmt_i] * hypothesis_scale - PEbyOpDet[pmt_i]) /
                 std::sqrt(lightHypothesis[pmt_i] * hypothesis_scale);

    hypothesis_integral += lightHypothesis[pmt_i] * hypothesis_scale;
  }

  *output << "\n\t TOTAL (hypothesis,flash)=(" << hypothesis_integral << "," << flash_integral
          << ")"
          << "  difference="
          << (hypothesis_integral - flash_integral) / std::sqrt(hypothesis_integral);

  *output << std::endl;
  *output << "End result=" << result << std::endl;
  *output << "----------------------------------------------" << std::endl;
}
