#include "larana/T0Finder/AssociationsTools/IHitParticleAssociations.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/MCParticle.h"

namespace t0 {
  ////////////////////////////////////////////////////////////////////////
  //
  // Class:       DirectHitParticleAssns
  // Module Type: art tool
  // File:        DirectHitParticleAssns.h
  //
  //              This provides MC truth information by using output
  //              reco Hit <--> MCParticle associations
  //
  // Configuration parameters:
  //
  // TruncMeanFraction     - the fraction of waveform bins to discard when
  //
  // Created by Tracy Usher (usher@slac.stanford.edu) on November 21, 2017
  //
  ////////////////////////////////////////////////////////////////////////

  class DirectHitParticleAssns : virtual public IHitParticleAssociations {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit DirectHitParticleAssns(fhicl::ParameterSet const& pset);

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const& pset) override;

    /**
     *  @brief This rebuilds the internal maps
     */
    void CreateHitParticleAssociations(art::Event&, HitParticleAssociations*) override;

  private:
    std::vector<art::InputTag> fHitModuleLabelVec;
    art::InputTag fMCParticleModuleLabel;

    struct TrackIDEinfo {
      float E;
      float NumElectrons;
    };
    std::unordered_map<int, TrackIDEinfo> fTrkIDECollector;
  };

  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  DirectHitParticleAssns::DirectHitParticleAssns(fhicl::ParameterSet const& pset)
  {
    reconfigure(pset);

    // Report.
    mf::LogInfo("DirectHitParticleAssns") << "Configured\n";
  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void DirectHitParticleAssns::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMCParticleModuleLabel = pset.get<art::InputTag>("MCParticleLabel");
    fHitModuleLabelVec = pset.get<std::vector<art::InputTag>>("HitModuleLabelVec");
  }

  //----------------------------------------------------------------------------
  /// Rebuild method -> rebuild the basic maps to get truth information
  ///
  /// Arguments:
  ///
  /// event - the art event used to extract all information
  ///
  void DirectHitParticleAssns::CreateHitParticleAssociations(art::Event& evt,
                                                             HitParticleAssociations* hitPartAssns)
  {
    // This function handles the "direct" creation of hit<-->MCParticle associations through use of the BackTracker
    //
    auto mcpartHandle = evt.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleModuleLabel);

    // Access art services...
    art::ServiceHandle<cheat::BackTrackerService const> btService;
    art::ServiceHandle<cheat::ParticleInventoryService const> piService;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Loop over input hit producer labels
    for (const auto& inputTag : fHitModuleLabelVec) {
      art::Handle<std::vector<recob::Hit>> hitListHandle;
      evt.getByLabel(inputTag, hitListHandle);

      if (!hitListHandle.isValid()) {
        mf::LogInfo("DirectHitParticleAssns")
          << "InputTag not associating to valid hit collection, tag: " << inputTag << "\n";
        continue;
      }

      anab::BackTrackerHitMatchingData bthmd;
      std::unordered_map<int, int>
        trkid_lookup; //indexed by geant4trkid, delivers MC particle location

      auto const& hitList(*hitListHandle);
      auto const& mcpartList(*mcpartHandle);

      for (size_t i_h = 0; i_h < hitList.size(); ++i_h) {
        art::Ptr<recob::Hit> hitPtr(hitListHandle, i_h);

        auto trkide_list = btService->HitToTrackIDEs(clockData, hitPtr);

        double maxe(-1.);
        double tote(0.);
        int maxtrkid(-1);
        double maxn(-1.);
        double totn(0.);
        int maxntrkid(-1);

        fTrkIDECollector.clear();

        //for(auto const& t : trkide_list){
        for (size_t i_t = 0; i_t < trkide_list.size(); ++i_t) {
          auto const& t(trkide_list[i_t]);
          fTrkIDECollector[t.trackID].E += t.energy;
          tote += t.energy;
          if (fTrkIDECollector[t.trackID].E > maxe) {
            maxe = fTrkIDECollector[t.trackID].E;
            maxtrkid = t.trackID;
          }
          fTrkIDECollector[t.trackID].NumElectrons += t.numElectrons;
          totn += t.numElectrons;
          if (fTrkIDECollector[t.trackID].NumElectrons > maxn) {
            maxn = fTrkIDECollector[t.trackID].NumElectrons;
            maxntrkid = t.trackID;
          }

          //if not found, find mc particle...
          if (trkid_lookup.find(t.trackID) == trkid_lookup.end()) {
            size_t i_p = 0;
            while (i_p < mcpartList.size()) {
              if (mcpartList[i_p].TrackId() == abs(t.trackID)) {
                trkid_lookup[t.trackID] = (int)i_p;
                break;
              }
              ++i_p;
            }
            if (i_p == mcpartList.size()) trkid_lookup[t.trackID] = -1;
          }
        }
        //end loop on TrackIDs

        //now find the mcparticle and loop back through ...
        for (auto const& t : fTrkIDECollector) {
          int mcpart_i = trkid_lookup[t.first];
          if (mcpart_i == -1) continue; //no mcparticle here
          art::Ptr<simb::MCParticle> mcpartPtr(mcpartHandle, mcpart_i);
          bthmd.ideFraction = t.second.E / tote;
          bthmd.isMaxIDE = (t.first == maxtrkid);
          bthmd.ideNFraction = t.second.NumElectrons / totn;
          bthmd.isMaxIDEN = (t.first == maxntrkid);
          bthmd.energy = t.second.E;
          bthmd.numElectrons = t.second.NumElectrons;
          hitPartAssns->addSingle(mcpartPtr, hitPtr, bthmd);
        }

      } //end loop on hits
    }   // end loop on producers

    return;
  }

  //----------------------------------------------------------------------------

  DEFINE_ART_CLASS_TOOL(DirectHitParticleAssns)
}
