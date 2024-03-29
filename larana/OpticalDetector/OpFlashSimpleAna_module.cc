////////////////////////////////////////////////////////////////////////
// Class:       OpFlashSimpleAna
// Module Type: analyzer
// File:        OpFlashSimpleAna_module.cc
//
// Generated at Tue Apr 21 10:58:24 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art_root_io/TFileService.h"

#include "OpFlashAnaAlg.h"

#include "TTree.h"

namespace opdet {
  class OpFlashSimpleAna;
}

class opdet::OpFlashSimpleAna : public art::EDAnalyzer {
public:
  explicit OpFlashSimpleAna(fhicl::ParameterSet const& p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpFlashSimpleAna(OpFlashSimpleAna const&) = delete;
  OpFlashSimpleAna(OpFlashSimpleAna&&) = delete;
  OpFlashSimpleAna& operator=(OpFlashSimpleAna const&) = delete;
  OpFlashSimpleAna& operator=(OpFlashSimpleAna&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;

private:
  // Declare member data here.
  std::string fOpFlashModuleLabel;
  std::string fOpHitModuleLabel;
  bool fMakeOpDetPEHist;

  OpFlashAnaAlg fAnaAlg;
};

opdet::OpFlashSimpleAna::OpFlashSimpleAna(fhicl::ParameterSet const& p) : EDAnalyzer(p) // ,
// More initializers here.
{
  fOpFlashModuleLabel = p.get<std::string>("OpFlashModuleLabel", "");
  fOpHitModuleLabel = p.get<std::string>("OpFlashModuleLabel", "");
  fMakeOpDetPEHist = p.get<bool>("MakeOpDetPEHist", true);
}

void opdet::OpFlashSimpleAna::analyze(art::Event const& e)
{
  if (fOpFlashModuleLabel.size() > 0) {
    art::Handle<std::vector<recob::OpFlash>> flashHandle;
    e.getByLabel(fOpFlashModuleLabel, flashHandle);
    std::vector<recob::OpFlash> const& flashVector(*flashHandle);
    fAnaAlg.FillOpFlashes(flashVector);
  }
  if (fOpHitModuleLabel.size() > 0) {
    art::Handle<std::vector<recob::OpHit>> hitHandle;
    e.getByLabel(fOpHitModuleLabel, hitHandle);
    std::vector<recob::OpHit> const& hitVector(*hitHandle);
    fAnaAlg.FillOpHits(hitVector);
  }
}

void opdet::OpFlashSimpleAna::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;
  if (fOpFlashModuleLabel.size() > 0)
    fAnaAlg.SetOpFlashTree(tfs->make<TTree>("OpFlashTree", "OpFlashSimpleAna: Flash Tree"),
                           fMakeOpDetPEHist);
  if (fOpHitModuleLabel.size() > 0)
    fAnaAlg.SetOpHitTree(tfs->make<TTree>("OpHitTree", "OpFlashSimpleAna: Hit Tree"));
}

DEFINE_ART_MODULE(opdet::OpFlashSimpleAna)
