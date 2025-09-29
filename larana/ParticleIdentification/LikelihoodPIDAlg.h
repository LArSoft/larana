////////////////////////////////////////////////////////////////////////
//
// A likelihood based particleID
//
// sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef LIKELIHOODPIDALG_H
#define LIKELIHOODPIDALG_H

#include <bitset>
#include <cmath>
#include <optional>
#include <string>

#include "fhiclcpp/fwd.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "larana/ParticleIdentification/PhysdEdx.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

class TProfile;

namespace anab {
  class Calorimetry;
  class ParticleID;
}

namespace pid {

  class LikelihoodPIDAlg {

  public:
    LikelihoodPIDAlg(fhicl::ParameterSet const& pset);

    /**
     * Helper function to go from geo::PlaneID to a bitset
     */
    std::bitset<8> GetBitset(geo::PlaneID planeID);

    anab::ParticleID DoParticleID(const std::vector<art::Ptr<anab::Calorimetry>>& calo);

  private:
    std::map<int, PhysdEdx*> map_PhysdEdx;
    trkf::TrackMomentumCalculator tmc;
    float fmaxrr;
  }; //
} // namespace
#endif // LIKELIHOODPIDALG_H
