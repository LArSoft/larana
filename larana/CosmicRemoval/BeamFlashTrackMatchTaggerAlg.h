#ifndef BEAMFLASHTRACKMATCHTAGGERALG_H
#define BEAMFLASHTRACKMATCHTAGGERALG_H
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

namespace detinfo {
  class LArProperties;
}
namespace phot {
  class PhotonVisibilityService;
}
namespace opdet {
  class OpDigiProperties;
}

#include "larcorealg/CoreUtils/ProviderPack.h"
#include "larcorealg/Geometry/fwd.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "fhiclcpp/fwd.h"

class TVector3;
class TH1F;
class TTree;

#include <iostream>
#include <string>
#include <vector>

namespace cosmic {
  class BeamFlashTrackMatchTaggerAlg;
}

class cosmic::BeamFlashTrackMatchTaggerAlg {
public:
  /// Pack of provider-interface supporting services we need
  using Providers_t = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;

  BeamFlashTrackMatchTaggerAlg(fhicl::ParameterSet const& p);

  //how to run the algorithm
  void RunCompatibilityCheck(std::vector<recob::OpFlash> const&,
                             std::vector<recob::Track> const&,
                             std::vector<anab::CosmicTag>&,
                             std::vector<size_t>&,
                             Providers_t,
                             phot::PhotonVisibilityService const&,
                             opdet::OpDigiProperties const&);

  void SetHypothesisComparisonTree(TTree*, TH1F*, TH1F*);

  void RunHypothesisComparison(unsigned int const,
                               unsigned int const,
                               std::vector<recob::OpFlash> const&,
                               std::vector<recob::Track> const&,
                               Providers_t,
                               phot::PhotonVisibilityService const&,
                               opdet::OpDigiProperties const&);

  void RunHypothesisComparison(unsigned int const,
                               unsigned int const,
                               std::vector<recob::OpFlash> const&,
                               std::vector<simb::MCParticle> const&,
                               Providers_t,
                               phot::PhotonVisibilityService const&,
                               opdet::OpDigiProperties const&);

private:
  const anab::CosmicTagID_t COSMIC_TYPE_FLASHMATCH;
  const anab::CosmicTagID_t COSMIC_TYPE_OUTSIDEDRIFT;
  const bool DEBUG_FLAG;

  float fMinTrackLength;
  float fMinOpHitPE;
  float fMIPdQdx;
  float fOpDetSaturation;
  float fSingleChannelCut;
  float fCumulativeChannelThreshold;
  unsigned int fCumulativeChannelCut;
  float fIntegralCut;

  bool fMakeOutsideDriftTags;
  bool fNormalizeHypothesisToFlash;

  TTree* cTree;

  typedef struct FlashComparisonProperties {
    unsigned int run;
    unsigned int event;

    unsigned int flash_index;
    float flash_totalPE;
    float flash_y;
    float flash_sigmay;
    float flash_z;
    float flash_sigmaz;
    unsigned int flash_nOpDet;

    unsigned int hyp_index;
    float hyp_totalPE;
    float hyp_y;
    float hyp_sigmay;
    float hyp_z;
    float hyp_sigmaz;

    float trk_startx;
    float trk_starty;
    float trk_startz;

    float trk_endx;
    float trk_endy;
    float trk_endz;

    float chi2;

    std::string leaf_structure;
    FlashComparisonProperties()
      : leaf_structure("run/i:event/i:flash_index/i:flash_totalPE/F:flash_y/F:flash_sigmay/"
                       "F:flash_z/F:flash_sigmaz/F:flash_nOpDet/i:hyp_index/i:hyp_totalPE/F:hyp_y/"
                       "F:hyp_sigmay/F:hyp_z/F:hyp_sigmaz/F:trk_startx/F:trk_starty/F:trk_startz/"
                       "F:trk_endx/F:trk_endy/F:trk_endz/F:chi2/F")
    {}

  } FlashComparisonProperties_t;

  FlashComparisonProperties_t cFlashComparison_p;
  std::vector<float> cOpDetVector_flash;
  std::vector<float> cOpDetVector_hyp;
  TH1F* cOpDetHist_flash;
  TH1F* cOpDetHist_hyp;

  typedef enum CompatibilityResultType {
    kCompatible = 0,
    kSingleChannelCut,
    kCumulativeChannelCut,
    kIntegralCut
  } CompatibilityResultType;

  //core functions
  std::vector<float> GetMIPHypotheses(recob::Track const& track,
                                      Providers_t providers,
                                      phot::PhotonVisibilityService const& pvs,
                                      opdet::OpDigiProperties const&,
                                      float XOffset = 0);

  std::vector<float> GetMIPHypotheses(simb::MCParticle const& particle,
                                      size_t start_i,
                                      size_t end_i,
                                      Providers_t providers,
                                      phot::PhotonVisibilityService const& pvs,
                                      opdet::OpDigiProperties const&,
                                      float XOffset = 0);

  void AddLightFromSegment(TVector3 const& pt1,
                           TVector3 const& pt2,
                           std::vector<float>& lightHypothesis,
                           float& totalHypothesisPE,
                           geo::GeometryCore const& geom,
                           phot::PhotonVisibilityService const& pvs,
                           float const& PromptMIPScintYield,
                           float XOffset);

  void NormalizeLightHypothesis(std::vector<float>& lightHypothesis,
                                float const& totalHypothesisPE,
                                geo::GeometryCore const& geom);

  CompatibilityResultType CheckCompatibility(std::vector<float> const& lightHypothesis,
                                             const recob::OpFlash* flashPointer,
                                             geo::GeometryCore const& geom,
                                             geo::WireReadoutGeom const& wireReadoutGeom);

  bool InDetector(TVector3 const&, geo::GeometryCore const&);
  bool InDriftWindow(double, double, geo::GeometryCore const&);

  void FillFlashProperties(std::vector<float> const& opdetVector,
                           float&,
                           float&,
                           float&,
                           float&,
                           float&,
                           geo::GeometryCore const& geom);

  float CalculateChi2(std::vector<float> const&, std::vector<float> const&);

  //debugging functions
  void PrintTrackProperties(recob::Track const&, std::ostream* output = &std::cout);
  void PrintFlashProperties(recob::OpFlash const&, std::ostream* output = &std::cout);
  void PrintHypothesisFlashComparison(std::vector<float> const&,
                                      const recob::OpFlash*,
                                      geo::GeometryCore const& geom,
                                      geo::WireReadoutGeom const& wireReadoutGeom,
                                      CompatibilityResultType,
                                      std::ostream* output = &std::cout);
};

#endif
