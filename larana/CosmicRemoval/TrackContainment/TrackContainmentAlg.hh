/**
 * \file TrackContainmentAlg.h
 *
 *
 * @author wketchum
*/

#ifndef TRK_TRACKCONTAINMENTALG_H
#define TRK_TRACKCONTAINMENTALG_H

#include "larcorealg/Geometry/fwd.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/Track.h"

#include "fhiclcpp/fwd.h"

class TTree;

#include <string>
#include <vector>

namespace trk {
  class TrackContainmentAlg;

  typedef struct TrackTree {

    TrackTree() {}

    TrackTree(recob::Track const& t)
      : start_x(t.Vertex().X())
      , start_y(t.Vertex().Y())
      , start_z(t.Vertex().Z())
      , start_theta(t.VertexDirection().Theta())
      , start_phi(t.VertexDirection().Phi())
      , end_x(t.End().X())
      , end_y(t.End().Y())
      , end_z(t.End().Z())
      , end_theta(t.EndDirection().Theta())
      , end_phi(t.EndDirection().Phi())
      , length(t.Length())
      , length_simple((t.End() - t.Vertex()).R())
    {}

    double start_x;
    double start_y;
    double start_z;
    double start_theta;
    double start_phi;
    double end_x;
    double end_y;
    double end_z;
    double end_theta;
    double end_phi;
    double length;
    double length_simple;

    std::string Leaflist()
    {
      return "start_x/D:start_y/D:start_z/D:start_theta/D:start_phi/D:end_x/D:end_y/D:end_z/"
             "D:end_theta/D:end_phi/D:length/D:length_simple/D";
    }

  } TrackTree_t;

}

class trk::TrackContainmentAlg {

public:
  /// Default constructor
  TrackContainmentAlg();

  void SetupOutputTree(TTree*);

  void Configure(fhicl::ParameterSet const&);

  void SetRunEvent(unsigned int const&, unsigned int const&);
  void ProcessTracks(std::vector<std::vector<recob::Track>> const&,
                     geo::GeometryCore const&,
                     geo::WireReadoutGeom const&);

  std::vector<std::vector<int>> const& GetTrackContainmentValues()
  {
    return fTrackContainmentLevel;
  }
  std::vector<std::vector<double>> const& GetTrackMinDistanceValues() { return fMinDistances; }
  std::vector<std::vector<anab::CosmicTag>> const& GetTrackCosmicTags();

  void setMakeCosmicTags(bool flag = true) { fMakeCosmicTags = flag; }
  void setFillOutputTree(bool flag = true) { fFillOutputTree = flag; }

private:
  double fZBuffer;
  double fYBuffer;
  double fXBuffer;
  double fIsolation;
  bool fMakeCosmicTags;
  bool fDebug;
  bool fFillOutputTree;

  TTree* fTrackTree;
  TrackTree_t fTrackTreeObj;
  unsigned int fRun;
  unsigned int fEvent;
  unsigned int fCollection;
  unsigned int fTrkID;
  double fDistance;
  int fContainment;

  std::vector<std::vector<int>> fTrackContainmentLevel;
  std::vector<std::vector<std::pair<int, int>>> fTrackContainmentIndices;
  std::vector<std::vector<double>> fMinDistances;
  std::vector<std::vector<anab::CosmicTag>> fCosmicTags;

  bool IsContained(recob::Track const&, geo::GeometryCore const&);
  anab::CosmicTagID_t GetCosmicTagID(recob::Track const&, geo::GeometryCore const&);

  double MinDistanceStartPt(recob::Track const&, recob::Track const&);
  double MinDistanceEndPt(recob::Track const&, recob::Track const&);
};

#endif
