////////////////////////////////////////////////////////////////////////
// Class:       CosmicPFParticleTagger
// Module Type: producer
// File:        CosmicPFParticleTagger_module.cc
//              This module checks timing and TPC volume boundaries as a
//              way to tag potential cosmic rays
//              This particular module uses PFParticles as input and handles
//              special cases associated with them.
//              This module started life as CosmicTrackTagger_module, written
//              by Sarah Lockwitz, and small alterations made to handle the
//              PFParticle input
//
// Generated at Wed Sep 17 19:17:00 2014 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicPFParticleTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <iterator>

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

namespace cosmic {
  class CosmicPFParticleTagger;
}

class cosmic::CosmicPFParticleTagger : public art::EDProducer {
public:
  explicit CosmicPFParticleTagger(fhicl::ParameterSet const& p);

  void produce(art::Event& e) override;

private:
  std::string fPFParticleModuleLabel;
  std::string fTrackModuleLabel;
  int fEndTickPadding;
  int fDetectorWidthTicks;
  int fMinTickDrift, fMaxTickDrift;
  int fMaxOutOfTime; ///< Max hits that can be out of time before rejecting
  float fTPCXBoundary, fTPCYBoundary, fTPCZBoundary;
  float fDetHalfHeight, fDetWidth, fDetLength;
};

cosmic::CosmicPFParticleTagger::CosmicPFParticleTagger(fhicl::ParameterSet const& p) : EDProducer{p}
{
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC({0, 0});
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);

  fDetHalfHeight = tpc.HalfHeight();
  fDetWidth = 2. * tpc.HalfWidth();
  fDetLength = tpc.Length();

  float fSamplingRate = sampling_rate(clock_data);

  fPFParticleModuleLabel = p.get<std::string>("PFParticleModuleLabel");
  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel", "track");
  fEndTickPadding = p.get<int>("EndTickPadding", 50); // Fudge the TPC edge in ticks...
  fMaxOutOfTime = p.get<int>("MaxOutOfTime", 4);

  fTPCXBoundary = p.get<float>("TPCXBoundary", 5);
  fTPCYBoundary = p.get<float>("TPCYBoundary", 5);
  fTPCZBoundary = p.get<float>("TPCZBoundary", 5);

  const double driftVelocity = detp.DriftVelocity(detp.Efield(), detp.Temperature()); // cm/us

  fDetectorWidthTicks =
    2 * tpc.HalfWidth() / (driftVelocity * fSamplingRate / 1000); // ~3200 for uB
  fMinTickDrift = clock_data.Time2Tick(clock_data.TriggerTime());
  fMaxTickDrift = fMinTickDrift + fDetectorWidthTicks + fEndTickPadding;

  produces<std::vector<anab::CosmicTag>>();
  produces<art::Assns<anab::CosmicTag, recob::Track>>();
  produces<art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

void cosmic::CosmicPFParticleTagger::produce(art::Event& evt)
{
  // Instatiate the output
  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagTrackVector(
    new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<anab::CosmicTag, recob::Track>> assnOutCosmicTagTrack(
    new art::Assns<anab::CosmicTag, recob::Track>);
  std::unique_ptr<art::Assns<recob::PFParticle, anab::CosmicTag>> assnOutCosmicTagPFParticle(
    new art::Assns<recob::PFParticle, anab::CosmicTag>);

  // Recover handle for PFParticles
  art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
  evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);

  if (!pfParticleHandle.isValid()) {
    evt.put(std::move(cosmicTagTrackVector));
    evt.put(std::move(assnOutCosmicTagTrack));
    return;
  }

  // Recover the handle for the tracks
  art::Handle<std::vector<recob::Track>> trackHandle;
  evt.getByLabel(fTrackModuleLabel, trackHandle);

  if (!trackHandle.isValid()) {
    evt.put(std::move(cosmicTagTrackVector));
    evt.put(std::move(assnOutCosmicTagTrack));
    return;
  }

  // Recover handle for track <--> PFParticle associations
  art::Handle<art::Assns<recob::PFParticle, recob::Track>> pfPartToTrackHandle;
  evt.getByLabel(fTrackModuleLabel, pfPartToTrackHandle);

  // Recover the list of associated tracks
  art::FindManyP<recob::Track> pfPartToTrackAssns(pfParticleHandle, evt, fTrackModuleLabel);

  // and the hits
  art::FindManyP<recob::Hit> hitsSpill(trackHandle, evt, fTrackModuleLabel);

  // The outer loop is going to be over PFParticles
  for (size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++) {
    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfPartIdx);

    // Recover the track vector
    std::vector<art::Ptr<recob::Track>> trackVec = pfPartToTrackAssns.at(pfPartIdx);

    // Is there a track associated to this PFParticle?
    if (trackVec.empty()) {
      // We need to make a null CosmicTag to store with this PFParticle to keep sequencing correct
      std::vector<float> tempPt1, tempPt2;
      tempPt1.push_back(-999);
      tempPt1.push_back(-999);
      tempPt1.push_back(-999);
      tempPt2.push_back(-999);
      tempPt2.push_back(-999);
      tempPt2.push_back(-999);
      cosmicTagTrackVector->emplace_back(tempPt1, tempPt2, 0., anab::CosmicTagID_t::kNotTagged);
      util::CreateAssn(evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
      continue;
    }

    // Start the tagging process...
    int isCosmic = 0;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;
    art::Ptr<recob::Track> track1 = trackVec.front();

    std::vector<art::Ptr<recob::Hit>> hitVec = hitsSpill.at(track1.key());

    // Recover track end points
    auto vertexPosition = track1->Vertex();
    auto vertexDirection = track1->VertexDirection();
    auto endPosition = track1->End();

    // In principle there is one track associated to a PFParticle... but with current
    // technology it can happen that a PFParticle is broken into multiple tracks. Our
    // aim here is to find the maximum extents of all the tracks which have been
    // associated to the single PFParticle
    if (trackVec.size() > 1) {
      for (size_t trackIdx = 1; trackIdx < trackVec.size(); trackIdx++) {
        art::Ptr<recob::Track> track(trackVec[trackIdx]);

        auto trackStart = track->Vertex();
        auto trackEnd = track->End();

        // Arc length possibilities for start of track
        double arcLStartToStart = (trackStart - vertexPosition).Dot(vertexDirection);
        double arcLStartToEnd = (trackEnd - vertexPosition).Dot(vertexDirection);

        if (arcLStartToStart < 0. || arcLStartToEnd < 0.) {
          if (arcLStartToStart < arcLStartToEnd)
            vertexPosition = trackStart;
          else
            vertexPosition = trackEnd;
        }

        // Arc length possibilities for end of track
        double arcLEndToStart = (trackStart - endPosition).Dot(vertexDirection);
        double arcLEndToEnd = (trackEnd - endPosition).Dot(vertexDirection);

        if (arcLEndToStart > 0. || arcLEndToEnd > 0.) {
          if (arcLEndToStart > arcLEndToEnd)
            endPosition = trackStart;
          else
            endPosition = trackEnd;
        }

        // add the hits from this track to the collection
        hitVec.insert(
          hitVec.end(), hitsSpill.at(track.key()).begin(), hitsSpill.at(track.key()).end());
      }
    }

    // "Track" end points in easily readable form
    float trackEndPt1_X = vertexPosition.X();
    float trackEndPt1_Y = vertexPosition.Y();
    float trackEndPt1_Z = vertexPosition.Z();
    float trackEndPt2_X = endPosition.X();
    float trackEndPt2_Y = endPosition.Y();
    float trackEndPt2_Z = endPosition.Z();

    /////////////////////////////////////
    // Check that all hits on particle are "in time"
    /////////////////////////////////////
    int nOutOfTime(0);

    for (unsigned int p = 0; p < hitVec.size(); p++) {
      int peakLessRms = hitVec[p]->PeakTimeMinusRMS();
      int peakPlusRms = hitVec[p]->PeakTimePlusRMS();

      if (peakLessRms < fMinTickDrift || peakPlusRms > fMaxTickDrift) {
        if (++nOutOfTime > fMaxOutOfTime) {
          isCosmic = 1;
          tag_id = anab::CosmicTagID_t::kOutsideDrift_Partial;
          break; // If one hit is out of time it must be a cosmic ray
        }
      }
    }

    /////////////////////////////////
    // Now check the TPC boundaries:
    /////////////////////////////////
    if (isCosmic == 0) {
      // In below we check entry and exit points. Note that a special case of a particle entering
      // and exiting the same surface is considered to be running parallel to the surface and NOT
      // entering and exiting.
      // Also, in what follows we make no assumptions on which end point is the "start" or
      // "end" of the track being considered.
      unsigned boundaryMask[] = {0, 0};

      // Check x extents - note that uboone coordinaes system has x=0 at edge
      // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
      // Also note that, in theory, any cosmic ray entering or exiting the X surfaces will have presumably
      // been removed already by the checking of "out of time" hits... but this will at least label
      // neutrino interaction tracks which exit through the X surfaces of the TPC
      if (fDetWidth - trackEndPt1_X < fTPCXBoundary)
        boundaryMask[0] = 0x1;
      else if (trackEndPt1_X < fTPCXBoundary)
        boundaryMask[0] = 0x2;

      if (fDetWidth - trackEndPt2_X < fTPCXBoundary)
        boundaryMask[1] = 0x1;
      else if (trackEndPt2_X < fTPCXBoundary)
        boundaryMask[1] = 0x2;

      // Check y extents (note coordinate system change)
      // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
      if (fDetHalfHeight - trackEndPt1_Y < fTPCYBoundary)
        boundaryMask[0] = 0x10;
      else if (fDetHalfHeight + trackEndPt1_Y < fTPCYBoundary)
        boundaryMask[0] = 0x20;

      if (fDetHalfHeight - trackEndPt2_Y < fTPCYBoundary)
        boundaryMask[1] = 0x10;
      else if (fDetHalfHeight + trackEndPt2_Y < fTPCYBoundary)
        boundaryMask[1] = 0x20;

      // Check z extents
      // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
      if (fDetLength - trackEndPt1_Z < fTPCZBoundary)
        boundaryMask[0] = 0x100;
      else if (trackEndPt1_Z < fTPCZBoundary)
        boundaryMask[0] = 0x200;

      if (fDetLength - trackEndPt2_Z < fTPCZBoundary)
        boundaryMask[1] = 0x100;
      else if (trackEndPt2_Z < fTPCZBoundary)
        boundaryMask[1] = 0x200;

      unsigned trackMask = boundaryMask[0] | boundaryMask[1];
      int nBitsSet(0);

      for (int idx = 0; idx < 12; idx++)
        if (trackMask & (0x1 << idx)) nBitsSet++;

      // This should check for the case of a track which is both entering and exiting
      // but we consider entering and exiting the z boundaries to be a special case (should it be?)
      if (nBitsSet > 1) {
        if ((trackMask & 0x300) != 0x300) {
          isCosmic = 2;
          if ((trackMask & 0x3) == 0x3)
            tag_id = anab::CosmicTagID_t::kGeometry_XX;
          else if ((trackMask & 0x30) == 0x30)
            tag_id = anab::CosmicTagID_t::kGeometry_YY;
          else if ((trackMask & 0x3) && (trackMask & 0x30))
            tag_id = anab::CosmicTagID_t::kGeometry_XY;
          else if ((trackMask & 0x3) && (trackMask & 0x300))
            tag_id = anab::CosmicTagID_t::kGeometry_XZ;
          else
            tag_id = anab::CosmicTagID_t::kGeometry_YZ;
        }
        // This is the special case of track which appears to enter/exit z boundaries
        else {
          isCosmic = 3;
          tag_id = anab::CosmicTagID_t::kGeometry_ZZ;
        }
      }
      // This looks for track which enters/exits a boundary but has other endpoint in TPC
      else if (nBitsSet > 0) {
        isCosmic = 4;
        if (trackMask & 0x3)
          tag_id = anab::CosmicTagID_t::kGeometry_X;
        else if (trackMask & 0x30)
          tag_id = anab::CosmicTagID_t::kGeometry_Y;
        else if (trackMask & 0x300)
          tag_id = anab::CosmicTagID_t::kGeometry_Z;
      }
    }

    std::vector<float> endPt1;
    std::vector<float> endPt2;
    endPt1.push_back(trackEndPt1_X);
    endPt1.push_back(trackEndPt1_Y);
    endPt1.push_back(trackEndPt1_Z);
    endPt2.push_back(trackEndPt2_X);
    endPt2.push_back(trackEndPt2_Y);
    endPt2.push_back(trackEndPt2_Z);

    float cosmicScore = isCosmic > 0 ? 1. : 0.;

    // Handle special cases
    if (isCosmic == 3)
      cosmicScore = 0.4; // Enter/Exit at opposite Z boundaries
    else if (isCosmic == 4)
      cosmicScore = 0.5; // Enter or Exit but not both

    // Loop through the tracks resulting from this PFParticle and mark them
    cosmicTagTrackVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id);

    util::CreateAssn(evt, *cosmicTagTrackVector, trackVec, *assnOutCosmicTagTrack);

    // Don't forget the association to the PFParticle
    util::CreateAssn(evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
  }

  evt.put(std::move(cosmicTagTrackVector));
  evt.put(std::move(assnOutCosmicTagTrack));
  evt.put(std::move(assnOutCosmicTagPFParticle));

} // end of produce

DEFINE_ART_MODULE(cosmic::CosmicPFParticleTagger)
