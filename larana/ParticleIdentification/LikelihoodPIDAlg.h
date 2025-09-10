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
#include <optional>
#include <string>
#include <cmath>

namespace fhicl {
  class ParameterSet;
}
#include "larana/ParticleIdentification/PhysdEdx.h"
#include "canvas/Persistency/Common/Ptr.h"

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
    std::map< int, PhysdEdx* > map_PhysdEdx;

    float M_mu = 105.65837; // MeV
    float M_pi = 139.570; // MeV, charged pion
    float M_pro = 938.272; // MeV
    
    trkf::TrackMomentumCalculator tmc;
    int fSkipNhits;

  }; //
} // namespace
#endif // LIKELIHOODPIDALG_H
