/////////////////////////////////////////////////////////////////////////////
/// Class:       MCParticleShowerMatching
/// Module Type: producer
/// File:        MCParticleShowerMatching_module.cc
///
/// Author:         Wesley Ketchum
/// E-mail address: wketchum@fnal.gov
///
/// This module uses existing hit<-->MCParticle assns to make assns of showers
/// to mcparticles. It also stores some idea of cleanliness of the matching.
/// Note: it's probably better to do fancier things with the hit<->particle
/// info, but this is a start of sorts.
///
/// Input: recob::Showers, recob::Hit collection, recob::Hit<--->simb::MCparticle assns
/// Output: recob::Shower/simb::MCParticle assns, with BackShowererMatchingData.
///
/////////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <iostream>
#include <iterator>
#include <memory>

// LArSoft
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace t0 {
  class MCParticleShowerMatching;
}

class t0::MCParticleShowerMatching : public art::EDProducer {
public:
  explicit MCParticleShowerMatching(fhicl::ParameterSet const& p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCParticleShowerMatching(MCParticleShowerMatching const&) = delete;
  MCParticleShowerMatching(MCParticleShowerMatching&&) = delete;
  MCParticleShowerMatching& operator=(MCParticleShowerMatching const&) = delete;
  MCParticleShowerMatching& operator=(MCParticleShowerMatching&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  art::InputTag fShowerModuleLabel;
  art::InputTag fShowerHitAssnLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fHitParticleAssnLabel;
};

t0::MCParticleShowerMatching::MCParticleShowerMatching(fhicl::ParameterSet const& p) : EDProducer{p}
{
  fShowerModuleLabel = p.get<art::InputTag>("ShowerModuleLabel");
  fShowerHitAssnLabel = p.get<art::InputTag>("ShowerHitAssnLabel", fShowerModuleLabel);
  fHitModuleLabel = p.get<art::InputTag>("HitModuleLabel");
  fHitParticleAssnLabel = p.get<art::InputTag>("HitParticleAssnLabel");

  produces<art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData>>();
}

void t0::MCParticleShowerMatching::produce(art::Event& evt)
{
  if (evt.isRealData()) return;

  //auto mcpartHandle = evt.getValidHandle< std::vector<simb::MCParticle> >("largeant");
  std::unique_ptr<art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData>>
    MCPartShowerassn(
      new art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData>);

  double maxe = -1;
  double tote = 0;
  // int    trkid = -1;
  //int    maxtrkid = -1;
  //double maxn = -1;
  //double totn = 0;
  //int maxntrkid = -1;

  anab::BackTrackerMatchingData btdata;
  std::unordered_map<int, double> trkide;

  art::Handle<std::vector<recob::Shower>> showerListHandle;
  evt.getByLabel(fShowerModuleLabel, showerListHandle);

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  evt.getByLabel(fHitModuleLabel, hitListHandle);

  if (!showerListHandle.isValid()) {
    std::cerr << "Shower handle is not valid!" << std::endl;
    return;
  }

  if (!hitListHandle.isValid()) {
    std::cerr << "Hit handle is not valid!" << std::endl;
    return;
  }

  auto const& showerList(*showerListHandle);
  art::FindManyP<recob::Hit> fmtht(showerListHandle, evt, fShowerHitAssnLabel);
  //auto const& mcpartList(*mcpartHandle);

  for (size_t i_t = 0; i_t < showerList.size(); ++i_t) {
    art::Ptr<recob::Shower> shwPtr(showerListHandle, i_t);
    trkide.clear();
    tote = 0;
    maxe = -1;
    art::Ptr<simb::MCParticle> maxp;

    std::vector<art::Ptr<recob::Hit>> allHits = fmtht.at(i_t);

    std::vector<anab::BackTrackerHitMatchingData const*> bthmd_vec;
    std::vector<art::Ptr<simb::MCParticle>> matchedParticlePtrs;

    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> particles_per_hit(
      hitListHandle, evt, fHitParticleAssnLabel);

    for (size_t i_h = 0; i_h < allHits.size(); ++i_h) {
      bthmd_vec.clear();
      matchedParticlePtrs.clear();
      particles_per_hit.get(allHits[i_h].key(), matchedParticlePtrs, bthmd_vec);

      for (size_t i_p = 0; i_p < matchedParticlePtrs.size(); ++i_p) {
        trkide[matchedParticlePtrs[i_p]->TrackId()] += bthmd_vec[i_p]->energy;
        tote += bthmd_vec[i_p]->energy;
        if (trkide[matchedParticlePtrs[i_p]->TrackId()] > maxe) {
          maxe = trkide[matchedParticlePtrs[i_p]->TrackId()];
          maxp = matchedParticlePtrs[i_p];
        }
      } //end loop over particles per hit

    } //end loop over hits

    btdata.cleanliness = maxe / tote;
    if (maxe > 0) MCPartShowerassn->addSingle(shwPtr, maxp, btdata);

  } //end loop over showers

  evt.put(std::move(MCPartShowerassn));
} // Produce

DEFINE_ART_MODULE(t0::MCParticleShowerMatching)
