////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Module Type: producer
// File:        CosmicTrackTagger_module.cc
//
// Generated at Mon Sep 24 18:21:00 2012 by Sarah Lockwitz using artmod
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicTrackTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <iostream>

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "TStopwatch.h"

namespace cosmic {
  class CosmicTrackTagger;
}

class cosmic::CosmicTrackTagger : public art::EDProducer {
public:
  explicit CosmicTrackTagger(fhicl::ParameterSet const& p);

  void produce(art::Event& e) override;

private:
  std::string fTrackModuleLabel;
  int fEndTickPadding;
  int fDetectorWidthTicks;
  float fTPCXBoundary, fTPCYBoundary, fTPCZBoundary;
  float fDetHalfHeight, fDetWidth, fDetLength;
  int fMinTickDrift, fMaxTickDrift;
};

cosmic::CosmicTrackTagger::CosmicTrackTagger(fhicl::ParameterSet const& p) : EDProducer{p}
{
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC({0, 0});

  fDetHalfHeight = tpc.HalfHeight();
  fDetWidth = 2. * tpc.HalfWidth();
  fDetLength = tpc.Length();

  float fSamplingRate = sampling_rate(clock_data);

  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel", "track");
  fEndTickPadding = p.get<int>("EndTickPadding", 50);

  fTPCXBoundary = p.get<float>("TPCXBoundary", 5);
  fTPCYBoundary = p.get<float>("TPCYBoundary", 5);
  fTPCZBoundary = p.get<float>("TPCZBoundary", 5);

  const double driftVelocity = detp.DriftVelocity(detp.Efield(), detp.Temperature()); // cm/us

  fDetectorWidthTicks =
    2 * tpc.HalfWidth() / (driftVelocity * fSamplingRate / 1000); // ~3200 for uB
  fMinTickDrift = clock_data.Time2Tick(clock_data.TriggerTime());
  fMaxTickDrift = fMinTickDrift + fDetectorWidthTicks + fEndTickPadding;

  produces<std::vector<anab::CosmicTag>>();
  produces<art::Assns<recob::Track, anab::CosmicTag>>();
}

void cosmic::CosmicTrackTagger::produce(art::Event& e)
{
  // Implementation of required member function here.

  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagTrackVector(
    new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<recob::Track, anab::CosmicTag>> assnOutCosmicTagTrack(
    new art::Assns<recob::Track, anab::CosmicTag>);

  TStopwatch ts;

  art::Handle<std::vector<recob::Track>> Trk_h;
  e.getByLabel(fTrackModuleLabel, Trk_h);
  std::vector<art::Ptr<recob::Track>> TrkVec;
  art::fill_ptr_vector(TrkVec, Trk_h);

  /////////////////////////////////
  // LOOPING OVER INSPILL TRACKS
  /////////////////////////////////

  art::FindManyP<recob::Hit> hitsSpill(Trk_h, e, fTrackModuleLabel);

  for (unsigned int iTrack = 0; iTrack < Trk_h->size(); iTrack++) {

    int isCosmic = 0;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;

    art::Ptr<recob::Track> tTrack = TrkVec.at(iTrack);
    std::vector<art::Ptr<recob::Hit>> HitVec = hitsSpill.at(iTrack);

    if (iTrack != tTrack.key()) { std::cout << "Mismatch in track index/key" << std::endl; }

    // A BETTER WAY OF FINDING END POINTS:
    auto tVector1 = tTrack->Vertex();
    auto tVector2 = tTrack->End();

    float trackEndPt1_X = tVector1.X();
    float trackEndPt1_Y = tVector1.Y();
    float trackEndPt1_Z = tVector1.Z();
    float trackEndPt2_X = tVector2.X();
    float trackEndPt2_Y = tVector2.Y();
    float trackEndPt2_Z = tVector2.Z();

    if (trackEndPt1_X != trackEndPt1_X || trackEndPt1_Y != trackEndPt1_Y ||
        trackEndPt1_Z != trackEndPt1_Z || trackEndPt2_X != trackEndPt2_X ||
        trackEndPt2_Y != trackEndPt2_Y || trackEndPt2_Z != trackEndPt2_Z) {
      std::cerr << "!!! FOUND A PROBLEM... the length is: " << tTrack->Length()
                << " np: " << tTrack->NumberTrajectoryPoints() << " id: " << tTrack->ID() << " "
                << tTrack << std::endl;
      std::vector<float> tempPt1, tempPt2;
      tempPt1.push_back(-999);
      tempPt1.push_back(-999);
      tempPt1.push_back(-999);
      tempPt2.push_back(-999);
      tempPt2.push_back(-999);
      tempPt2.push_back(-999);
      cosmicTagTrackVector->emplace_back(tempPt1, tempPt2, -999, tag_id);
      util::CreateAssn(e, *cosmicTagTrackVector, tTrack, *assnOutCosmicTagTrack);
      continue; // I don't want to deal with these "tracks"
    }

    /////////////////////////////////////
    // Getting first and last ticks
    /////////////////////////////////////
    float tick1 = 9999;
    float tick2 = -9999;
    bool dumpMe(false);

    for (unsigned int p = 0; p < HitVec.size(); p++) {
      if (dumpMe) {
        std::cout << "###>> Hit key: " << HitVec[p].key()
                  << ", peak - RMS: " << HitVec[p]->PeakTimeMinusRMS()
                  << ", peak + RMS: " << HitVec[p]->PeakTimePlusRMS() << std::endl;
      }
      if (HitVec[p]->PeakTimeMinusRMS() < tick1) tick1 = HitVec[p]->PeakTimeMinusRMS();
      if (HitVec[p]->PeakTimePlusRMS() > tick2) tick2 = HitVec[p]->PeakTimePlusRMS();
    }

    /////////////////////////////////////////////////////////
    // Are any of the ticks outside of the ReadOutWindow ?
    /////////////////////////////////////////////////////////
    if (tick1 < fMinTickDrift || tick2 > fMaxTickDrift) {
      isCosmic = 1;
      tag_id = anab::CosmicTagID_t::kOutsideDrift_Partial;
    }

    /////////////////////////////////
    // Now check Y & Z boundaries:
    /////////////////////////////////
    int nBdY = 0, nBdZ = 0;
    if (isCosmic == 0) {

      // Checking lower side of TPC
      if (fabs(fDetHalfHeight + trackEndPt1_Y) < fTPCYBoundary ||
          fabs(fDetHalfHeight + trackEndPt2_Y) < fTPCYBoundary || trackEndPt1_Y < -fDetHalfHeight ||
          trackEndPt2_Y < -fDetHalfHeight)
        nBdY++;

      // Checking upper side of TPC
      if (fabs(fDetHalfHeight - trackEndPt1_Y) < fTPCYBoundary ||
          fabs(fDetHalfHeight - trackEndPt2_Y) < fTPCYBoundary || trackEndPt1_Y > fDetHalfHeight ||
          trackEndPt2_Y > fDetHalfHeight)
        nBdY++;

      if (fabs(trackEndPt1_Z - fDetLength) < fTPCZBoundary ||
          fabs(trackEndPt2_Z - fDetLength) < fTPCZBoundary)
        nBdZ++;
      if (fabs(trackEndPt1_Z) < fTPCZBoundary || fabs(trackEndPt2_Z) < fTPCZBoundary) nBdZ++;
      if ((nBdY + nBdZ) > 1) {
        isCosmic = 2;
        if (nBdY > 1)
          tag_id = anab::CosmicTagID_t::kGeometry_YY;
        else if (nBdZ > 1)
          tag_id = anab::CosmicTagID_t::kGeometry_ZZ;
        else
          tag_id = anab::CosmicTagID_t::kGeometry_YZ;
      }
      else if ((nBdY + nBdZ) == 1) {
        isCosmic = 3;
        if (nBdY == 1)
          tag_id = anab::CosmicTagID_t::kGeometry_Y;
        else if (nBdZ == 1)
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

    float cosmicScore = isCosmic > 0 ? 1 : 0;
    if (isCosmic == 3) cosmicScore = 0.5;

    ///////////////////////////////////////////////////////
    // Doing a very basic check on X boundaries
    // this gets the types of tracks that go through both X boundaries of the detector
    if (fabs(trackEndPt1_X - trackEndPt2_X) > fDetWidth - fTPCXBoundary) {
      cosmicScore = 1;
      isCosmic = 4;
      tag_id = anab::CosmicTagID_t::kGeometry_XX;
    }

    cosmicTagTrackVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id);

    util::CreateAssn(e, *cosmicTagTrackVector, tTrack, *assnOutCosmicTagTrack);
  }
  // END OF LOOPING OVER INSPILL TRACKS

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////TAGGING DELTA RAYS (and other stub) ASSOCIATED TO A ALREADY TAGGED COSMIC TRACK//////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  float dE = 0, dS = 0, temp = 0, IScore = 0;
  unsigned int IndexE = 0, iTrk1 = 0, iTrk = 0;
  anab::CosmicTagID_t IType = anab::CosmicTagID_t::kNotTagged;

  for (iTrk = 0; iTrk < Trk_h->size(); iTrk++) {
    art::Ptr<recob::Track> tTrk = TrkVec.at(iTrk);
    if ((*cosmicTagTrackVector)[iTrk].CosmicScore() == 0) {
      auto tStart = tTrk->Vertex();
      auto tEnd = tTrk->End();
      unsigned int l = 0;
      for (iTrk1 = 0; iTrk1 < Trk_h->size(); iTrk1++) {
        art::Ptr<recob::Track> tTrk1 = TrkVec.at(iTrk1);
        float getScore = (*cosmicTagTrackVector)[iTrk1].CosmicScore();
        if (getScore == 1 || getScore == 0.5) {
          anab::CosmicTagID_t getType = (*cosmicTagTrackVector)[iTrk1].CosmicType();
          auto tStart1 = tTrk1->Vertex();
          auto tEnd1 = tTrk1->End();
          auto NumE = (tEnd - tStart1).Cross(tEnd - tEnd1);
          auto DenE = tEnd1 - tStart1;
          dE = NumE.R() / DenE.R();
          if (l == 0) {
            temp = dE;
            IndexE = iTrk1;
            IScore = getScore;
            IType = getType;
          }
          if (dE < temp) {
            temp = dE;
            IndexE = iTrk1;
            IScore = getScore;
            IType = getType;
          }
          l++;
        }
      } //End Trk1 loop
      art::Ptr<recob::Track> tTrkI = TrkVec.at(IndexE);
      auto tStartI = tTrkI->Vertex();
      auto tEndI = tTrkI->End();
      auto NumS = (tStart - tStartI).Cross(tStart - tEndI);
      auto DenS = tEndI - tStartI;
      dS = NumS.R() / DenS.R();
      if (((dS < 5 && temp < 5) || (dS < temp && dS < 5)) && (tTrk->Length() < 60)) {
        (*cosmicTagTrackVector)[iTrk].CosmicScore() = IScore - 0.05;
        (*cosmicTagTrackVector)[iTrk].CosmicType() = IType;
      }
    } // end cosmicScore==0 loop
  }   // end iTrk loop

  e.put(std::move(cosmicTagTrackVector));
  e.put(std::move(assnOutCosmicTagTrack));

  TrkVec.clear();

} // end of produce

DEFINE_ART_MODULE(cosmic::CosmicTrackTagger)
