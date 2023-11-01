#include "TrackContainmentAlg.hh"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

#include <iostream>

trk::TrackContainmentAlg::TrackContainmentAlg() {}

void trk::TrackContainmentAlg::SetupOutputTree(TTree* tfs_tree_trk)
{
  fTrackTree = tfs_tree_trk;
  fTrackTree->SetObject(fTrackTree->GetName(), "Track Tree");
  fTrackTree->Branch("run", &fRun);
  fTrackTree->Branch("event", &fEvent);
  fTrackTree->Branch("trk", &fTrackTreeObj, fTrackTreeObj.Leaflist().c_str());
  fTrackTree->Branch("trk_col", &fCollection);
  fTrackTree->Branch("trk_id", &fTrkID);
  fTrackTree->Branch("trk_mindist", &fDistance);
  fTrackTree->Branch("trk_containment", &fContainment);
}

void trk::TrackContainmentAlg::Configure(fhicl::ParameterSet const& p)
{
  fZBuffer = p.get<double>("ZBuffer");
  fYBuffer = p.get<double>("YBuffer");
  fXBuffer = p.get<double>("XBuffer");
  fIsolation = p.get<double>("Isolation");

  fDebug = p.get<bool>("Debug", false);
  fMakeCosmicTags = p.get<bool>("MakeCosmicTags", true);
  fFillOutputTree = p.get<bool>("FillOutputTree", true);
}

bool trk::TrackContainmentAlg::IsContained(recob::Track const& track, geo::GeometryCore const& geo)
{
  auto const& tpc = geo.TPC({0, 0});
  if (track.Vertex().Z() < (0 + fZBuffer) || track.Vertex().Z() > (tpc.Length() - fZBuffer))
    return false;
  if (track.End().Z() < (0 + fZBuffer) || track.End().Z() > (tpc.Length() - fZBuffer)) return false;
  if (track.Vertex().Y() < (-1 * tpc.HalfHeight() + fYBuffer) ||
      track.Vertex().Y() > (tpc.HalfHeight() - fYBuffer))
    return false;
  if (track.End().Y() < (-1 * tpc.HalfHeight() + fYBuffer) ||
      track.End().Y() > (tpc.HalfHeight() - fYBuffer))
    return false;
  if (track.Vertex().X() < (0 + fXBuffer) || track.Vertex().X() > (2 * tpc.HalfWidth() - fXBuffer))
    return false;
  if (track.End().X() < (0 + fXBuffer) || track.End().X() > (2 * tpc.HalfWidth() - fXBuffer))
    return false;

  return true;
}

anab::CosmicTagID_t trk::TrackContainmentAlg::GetCosmicTagID(recob::Track const& track,
                                                             geo::GeometryCore const& geo)
{
  auto const& tpc = geo.TPC({0, 0});

  auto id = anab::CosmicTagID_t::kNotTagged;

  if (track.Vertex().Z() < (0 + fZBuffer) || track.Vertex().Z() > (tpc.Length() - fZBuffer))
    id = anab::CosmicTagID_t::kGeometry_Z;
  else if (track.Vertex().Y() < (-1 * tpc.HalfHeight() + fYBuffer) ||
           track.Vertex().Y() > (tpc.HalfHeight() - fYBuffer))
    id = anab::CosmicTagID_t::kGeometry_Y;
  else if ((track.Vertex().X() > 0 && track.Vertex().X() < (0 + fXBuffer)) ||
           (track.Vertex().X() < 2 * tpc.HalfWidth() &&
            track.Vertex().X() > (2 * tpc.HalfWidth() - fXBuffer)))
    id = anab::CosmicTagID_t::kGeometry_X;

  if (track.End().Z() < (0 + fZBuffer) || track.End().Z() > (tpc.Length() - fZBuffer)) {
    if (id == anab::CosmicTagID_t::kNotTagged)
      id = anab::CosmicTagID_t::kGeometry_Z;
    else if (id == anab::CosmicTagID_t::kGeometry_Z)
      id = anab::CosmicTagID_t::kGeometry_ZZ;
    else if (id == anab::CosmicTagID_t::kGeometry_Y)
      id = anab::CosmicTagID_t::kGeometry_YZ;
    else if (id == anab::CosmicTagID_t::kGeometry_X)
      id = anab::CosmicTagID_t::kGeometry_XZ;
  }
  else if (track.End().Y() < (-1 * tpc.HalfHeight() + fYBuffer) ||
           track.End().Y() > (tpc.HalfHeight() - fYBuffer)) {
    if (id == anab::CosmicTagID_t::kNotTagged)
      id = anab::CosmicTagID_t::kGeometry_Y;
    else if (id == anab::CosmicTagID_t::kGeometry_Z)
      id = anab::CosmicTagID_t::kGeometry_YZ;
    else if (id == anab::CosmicTagID_t::kGeometry_Y)
      id = anab::CosmicTagID_t::kGeometry_YY;
    else if (id == anab::CosmicTagID_t::kGeometry_X)
      id = anab::CosmicTagID_t::kGeometry_XY;
  }
  else if ((track.End().X() > 0 && track.End().X() < (0 + fXBuffer)) ||
           (track.End().X() < 2 * tpc.HalfWidth() &&
            track.End().X() > (2 * tpc.HalfWidth() - fXBuffer))) {
    if (id == anab::CosmicTagID_t::kNotTagged)
      id = anab::CosmicTagID_t::kGeometry_X;
    else if (id == anab::CosmicTagID_t::kGeometry_Z)
      id = anab::CosmicTagID_t::kGeometry_XZ;
    else if (id == anab::CosmicTagID_t::kGeometry_Y)
      id = anab::CosmicTagID_t::kGeometry_XY;
    else if (id == anab::CosmicTagID_t::kGeometry_X)
      id = anab::CosmicTagID_t::kGeometry_XX;
  }

  if (track.Vertex().X() < 0 || track.Vertex().X() > 2 * tpc.HalfWidth())
    id = anab::CosmicTagID_t::kOutsideDrift_Partial;
  if (track.End().X() < 0 || track.End().X() > 2 * tpc.HalfWidth()) {
    if (id == anab::CosmicTagID_t::kOutsideDrift_Partial)
      id = anab::CosmicTagID_t::kOutsideDrift_Complete;
    else
      id = anab::CosmicTagID_t::kOutsideDrift_Partial;
  }

  return id;
}

double trk::TrackContainmentAlg::MinDistanceStartPt(recob::Track const& t_probe,
                                                    recob::Track const& t_ref)
{
  double min_distance = 9e12;
  for (size_t i_p = 0; i_p < t_ref.NumberTrajectoryPoints(); ++i_p) {
    double const distance = (t_probe.Vertex() - t_ref.LocationAtPoint(i_p)).R();
    if (distance < min_distance) min_distance = distance;
  }
  return std::sqrt(min_distance);
}

double trk::TrackContainmentAlg::MinDistanceEndPt(recob::Track const& t_probe,
                                                  recob::Track const& t_ref)
{
  double min_distance = 9e12;
  for (size_t i_p = 0; i_p < t_ref.NumberTrajectoryPoints(); ++i_p) {
    double const distance = (t_probe.End() - t_ref.LocationAtPoint(i_p)).R();
    if (distance < min_distance) min_distance = distance;
  }
  return std::sqrt(min_distance);
}

void trk::TrackContainmentAlg::SetRunEvent(unsigned int const& run, unsigned int const& event)
{
  fRun = run;
  fEvent = event;
}

void trk::TrackContainmentAlg::ProcessTracks(
  std::vector<std::vector<recob::Track>> const& tracksVec,
  geo::GeometryCore const& geo,
  geo::WireReadoutGeom const& wireReadoutGeom)
{
  auto const& tpc = geo.TPC({0, 0});

  if (fDebug) {
    std::cout << "Geometry:" << std::endl;
    std::cout << "\t" << tpc.HalfWidth() << " " << tpc.HalfHeight() << " " << tpc.Length()
              << std::endl;
    std::cout << "\t z:(" << fZBuffer << "," << tpc.Length() - fZBuffer << ")"
              << "\t y:(" << -1. * tpc.HalfHeight() + fYBuffer << "," << tpc.HalfHeight() - fYBuffer
              << ")"
              << "\t x:(" << 0 + fXBuffer << "," << 2 * tpc.HalfWidth() - fXBuffer << ")"
              << std::endl;
  }

  int containment_level = 0;
  bool track_linked = false;

  fTrackContainmentLevel.clear();
  fTrackContainmentLevel.resize(tracksVec.size());
  fMinDistances.clear();
  fMinDistances.resize(tracksVec.size());

  fTrackContainmentIndices.clear();
  fTrackContainmentIndices.push_back(std::vector<std::pair<int, int>>());

  fCosmicTags.clear();
  fCosmicTags.resize(tracksVec.size());

  //first, loop through tracks and see what's not contained

  for (size_t i_tc = 0; i_tc < tracksVec.size(); ++i_tc) {
    fTrackContainmentLevel[i_tc].resize(tracksVec[i_tc].size(), -1);
    fMinDistances[i_tc].resize(tracksVec[i_tc].size(), 9e12);
    fCosmicTags[i_tc].resize(tracksVec[i_tc].size(), anab::CosmicTag(-1));
    for (size_t i_t = 0; i_t < tracksVec[i_tc].size(); ++i_t) {

      if (!IsContained(tracksVec[i_tc][i_t], geo)) {
        if (!track_linked) track_linked = true;
        fTrackContainmentLevel[i_tc][i_t] = 0;
        fTrackContainmentIndices.back().emplace_back(i_tc, i_t);
        if (fDebug) {
          std::cout << "\tTrack (" << i_tc << "," << i_t << ")"
                    << " " << containment_level << std::endl;
        }

      } //end if contained
    }   //end loop over tracks

  } //end loop over track collections

  //now, while we are still linking tracks, loop over all tracks and note anything
  //close to an uncontained (or linked) track
  while (track_linked) {

    track_linked = false;
    ++containment_level;
    fTrackContainmentIndices.push_back(std::vector<std::pair<int, int>>());

    for (size_t i_tc = 0; i_tc < tracksVec.size(); ++i_tc) {
      for (size_t i_t = 0; i_t < tracksVec[i_tc].size(); ++i_t) {
        if (fTrackContainmentLevel[i_tc][i_t] >= 0)
          continue;
        else {
          for (auto const& i_tr : fTrackContainmentIndices[containment_level - 1]) {

            if (MinDistanceStartPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]) <
                fMinDistances[i_tc][i_t])
              fMinDistances[i_tc][i_t] =
                MinDistanceStartPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]);
            if (MinDistanceEndPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]) <
                fMinDistances[i_tc][i_t])
              fMinDistances[i_tc][i_t] =
                MinDistanceEndPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]);

            if (MinDistanceStartPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]) <
                  fIsolation ||
                MinDistanceEndPt(tracksVec[i_tc][i_t], tracksVec[i_tr.first][i_tr.second]) <
                  fIsolation) {
              if (!track_linked) track_linked = true;
              fTrackContainmentLevel[i_tc][i_t] = containment_level;
              fTrackContainmentIndices.back().emplace_back(i_tc, i_t);

              if (fDebug) {
                std::cout << "\tTrackPair (" << i_tc << "," << i_t << ") and (" << i_tr.first << ","
                          << i_tr.second << ")"
                          << " " << containment_level << std::endl;
              }

            } //end if track not isolated

          } //end loop over existing uncontained/linked tracks

        } //end if track not already linked

      } //end loops over tracks
    }   //end loop over track collections
  }     //end while linking tracks

  if (fDebug) std::cout << "All done! Now let's make the tree and tags!" << std::endl;

  //now we're going to will the tree and create tags if requested
  for (size_t i_tc = 0; i_tc < tracksVec.size(); ++i_tc) {
    for (size_t i_t = 0; i_t < tracksVec[i_tc].size(); ++i_t) {

      //fill ROOT Tree
      if (fFillOutputTree) {
        fTrackTreeObj = TrackTree_t(tracksVec[i_tc][i_t]);
        fDistance = fMinDistances[i_tc][i_t];
        fCollection = i_tc;
        fTrkID = i_t;
        fContainment = fTrackContainmentLevel[i_tc][i_t];
        fTrackTree->Fill();
      }

      if (fMakeCosmicTags) {

        //default (if track looks contained and isolated)
        float score = 0;
        auto id = anab::CosmicTagID_t::kNotTagged;

        //overwrite if track is not contained or not isolated
        if (fTrackContainmentLevel[i_tc][i_t] >= 0) {
          score = 1. / (1. + (float)fTrackContainmentLevel[i_tc][i_t]);
          if (fTrackContainmentLevel[i_tc][i_t] == 0)
            id = GetCosmicTagID(tracksVec[i_tc][i_t], geo);
          else
            id = anab::CosmicTagID_t::kNotIsolated;
        }

        fCosmicTags[i_tc][i_t] =
          anab::CosmicTag(std::vector<float>{(float)tracksVec[i_tc][i_t].Vertex().X(),
                                             (float)tracksVec[i_tc][i_t].Vertex().Y(),
                                             (float)tracksVec[i_tc][i_t].Vertex().Z()},
                          std::vector<float>{(float)tracksVec[i_tc][i_t].End().X(),
                                             (float)tracksVec[i_tc][i_t].End().Y(),
                                             (float)tracksVec[i_tc][i_t].End().Z()},
                          score,
                          id);
      } //end cosmic tag making

      //some debug work
      if (fTrackContainmentLevel[i_tc][i_t] < 0 && fDebug) {
        std::cout << "Track (" << i_tc << "," << i_t << ")"
                  << " " << fTrackContainmentLevel[i_tc][i_t] << " " << fMinDistances[i_tc][i_t]
                  << std::endl;
        auto const& vertex = tracksVec[i_tc][i_t].Vertex();
        std::cout << "\tS_(X,Y,Z) = (" << vertex.X() << "," << vertex.Y() << "," << vertex.Z()
                  << ")\n";
        std::cout << "\tNearest wire ..." << std::endl;
        for (unsigned int i_p = 0; i_p < wireReadoutGeom.Nplanes(); ++i_p) {
          auto const& plane = wireReadoutGeom.Plane({0, 0, i_p});
          std::cout << "\t\tPlane " << i_p << " " << plane.NearestWireID(vertex).Wire << std::endl;
        }

        auto const& end = tracksVec[i_tc][i_t].End();
        std::cout << "\tE_(X,Y,Z) = (" << end.X() << "," << end.Y() << "," << end.Z() << ")\n";
        std::cout << "\tNearest wire ..." << std::endl;
        for (unsigned int i_p = 0; i_p < wireReadoutGeom.Nplanes(); ++i_p) {
          auto const& plane = wireReadoutGeom.Plane({0, 0, i_p});
          std::cout << "\t\tPlane " << i_p << " " << plane.NearestWireID(end).Wire << std::endl;
        }
        std::cout << "\tLength=" << tracksVec[i_tc][i_t].Length() << std::endl;
        std::cout << "\tSimple_length=" << (end - vertex).R() << std::endl;
      } //end debug statements if track contained

    } //end loops over tracks
  }   //end loop over track collections

} //end ProcessTracks

std::vector<std::vector<anab::CosmicTag>> const& trk::TrackContainmentAlg::GetTrackCosmicTags()
{
  if (!fMakeCosmicTags)
    throw cet::exception("TrackContainmentAlg::GetTrackCosmicTags")
      << "Cosmic tags not created. Set MakeCosmicTags to true in fcl paramters.";
  return fCosmicTags;
}
