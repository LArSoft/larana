////////////////////////////////////////////////////////////////////////////
//
// \brief A likelihood based particle identification method using calorimetry information
//
// \author sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#include "larana/ParticleIdentification/LikelihoodPIDAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

namespace pid {
  class LikelihoodParticleID;
}

class pid::LikelihoodParticleID : public art::EDProducer {
public:
  explicit LikelihoodParticleID(fhicl::ParameterSet const& p);

  virtual void produce(art::Event& e);

private:
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;

  LikelihoodPIDAlg fLikelihoodAlg;
};

pid::LikelihoodParticleID::LikelihoodParticleID(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fTrackModuleLabel(p.get<std::string>("TrackModuleLabel"))
  , fCalorimetryModuleLabel(p.get<std::string>("CalorimetryModuleLabel"))
  , fLikelihoodAlg(p.get<fhicl::ParameterSet>("LikelihoodPIDAlg"))
{
  produces<std::vector<anab::ParticleID>>();
  produces<art::Assns<recob::Track, anab::ParticleID>>();
}

void pid::LikelihoodParticleID::produce(art::Event& evt)
{
  art::Handle<std::vector<recob::Track>> trackListHandle;
  evt.getByLabel(fTrackModuleLabel, trackListHandle);

  std::vector<art::Ptr<recob::Track>> tracklist;
  art::fill_ptr_vector(tracklist, trackListHandle);

  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  //if (!fmcal.isValid()) return;

  std::unique_ptr<std::vector<anab::ParticleID>> particleidcol(new std::vector<anab::ParticleID>);
  std::unique_ptr<art::Assns<recob::Track, anab::ParticleID>> assn(
    new art::Assns<recob::Track, anab::ParticleID>);

  if (fmcal.isValid()) {
    std::vector<art::Ptr<anab::Calorimetry>> calovec(1, art::Ptr<anab::Calorimetry>());
    for (size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter) {
      for (size_t i = 0; i < fmcal.at(trkIter).size(); ++i) {
        calovec[0] = fmcal.at(trkIter)[i];
        anab::ParticleID pidout = fLikelihoodAlg.DoParticleID(calovec);
        particleidcol->push_back(pidout);
        util::CreateAssn(evt, *particleidcol, tracklist[trkIter], *assn);
      }
    }
  }
  evt.put(std::move(particleidcol));
  evt.put(std::move(assn));

  return;
}

DEFINE_ART_MODULE(pid::LikelihoodParticleID)
