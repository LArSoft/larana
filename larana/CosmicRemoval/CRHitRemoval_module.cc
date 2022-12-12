////////////////////////////////////////////////////////////////////////
//
// Class:       CRHitRemoval
// Module Type: producer
// File:        CRHitRemoval_module.cc
//
// This module produces RecoBase/Hit objects after removing those
// deemed to be due to CR muons.
//
// Configuration parameters:
//
// CosmicProducerLabels    - a list of cosmic ray producers which should be or'ed
// HitProducerLabel        - the producer of the recob::Hit objects
// PFParticleProducerLabel - the producer of the recob::PFParticles to consider
// TrackProducerLabel      - the producer of the recob::Track objects
// CosmicTagThresholds     - a vector of thresholds to apply to label as cosmic
// EndTickPadding          - # ticks to "pad" the end tick to account for possible
//                           uncertainty in drift velocity
//
// Created by Tracy Usher (usher@slac.stanford.edu) on September 18, 2014
//
////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

class CRHitRemoval : public art::EDProducer {
public:
  // Copnstructors, destructor.
  explicit CRHitRemoval(fhicl::ParameterSet const& pset);

  // Overrides.
  virtual void produce(art::Event& e);
  virtual void beginJob();
  virtual void endJob();

private:
  // define vector for hits to make sure of uniform use
  using HitPtrVector = std::vector<art::Ptr<recob::Hit>>;

  // Methods
  void collectPFParticleHits(const recob::PFParticle* pfParticle,
                             const art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                             const art::FindManyP<recob::Cluster>& partToClusAssns,
                             const art::FindManyP<recob::Hit>& clusToHitAssns,
                             HitPtrVector& hitVec);

  void copyAllHits(std::vector<art::Ptr<recob::Hit>>&,
                   art::FindOneP<recob::Wire>&,
                   recob::HitCollectionCreator&);

  void copyInTimeHits(std::vector<art::Ptr<recob::Hit>>&,
                      art::FindOneP<recob::Wire>&,
                      recob::HitCollectionCreator&);

  void FilterHits(HitPtrVector& hits, HitPtrVector& used_hits);

  // Fcl parameters.
  std::vector<std::string> fCosmicProducerLabels; ///< List of cosmic tagger producers
  std::string fHitProducerLabel;                  ///< The full collection of hits
  std::string fPFParticleProducerLabel;           ///< PFParticle producer
  std::vector<std::string> fTrackProducerLabels;  ///< Track producer
  std::vector<std::string> fAssnProducerLabels;   ///< Track to PFParticle assns producer

  std::vector<double> fCosmicTagThresholds; ///< Thresholds for tagging

  int fEndTickPadding; ///< Padding the end tick

  int fDetectorWidthTicks; ///< Effective drift time in ticks
  int fMinTickDrift;       ///< Starting tick
  int fMaxTickDrift;       ///< Ending tick
  int fMaxOutOfTime;       ///< Max hits that can be out of time before rejecting

  // Statistics.
  int fNumEvent;     ///< Number of events seen.
  int fNumCRRejects; ///< Number of tracks produced.
};

DEFINE_ART_MODULE(CRHitRemoval)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CRHitRemoval::CRHitRemoval(fhicl::ParameterSet const& pset)
  : EDProducer{pset}, fNumEvent(0), fNumCRRejects(0)
{
  fCosmicProducerLabels = pset.get<std::vector<std::string>>("CosmicProducerLabels");
  fHitProducerLabel = pset.get<std::string>("HitProducerLabel");
  fPFParticleProducerLabel = pset.get<std::string>("PFParticleProducerLabel");
  fTrackProducerLabels = pset.get<std::vector<std::string>>("TrackProducerLabels");
  fAssnProducerLabels = pset.get<std::vector<std::string>>("AssnProducerLabels");
  fCosmicTagThresholds = pset.get<std::vector<double>>("CosmicTagThresholds");
  fEndTickPadding = pset.get<int>("EndTickPadding", 50);
  fMaxOutOfTime = pset.get<int>("MaxOutOfTime", 4);

  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
  recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);

  // Report.
  mf::LogInfo("CRHitRemoval") << "CRHitRemoval configured\n";
}

//----------------------------------------------------------------------------
/// Begin job method.
void CRHitRemoval::beginJob()
{
  auto const& tpc = lar::providerFrom<geo::Geometry>()->TPC({0, 0});
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);

  float const samplingRate = sampling_rate(clock_data);
  float const driftVelocity = detp.DriftVelocity(detp.Efield(), detp.Temperature()); // cm/us

  fDetectorWidthTicks = 2 * tpc.HalfWidth() / (driftVelocity * samplingRate / 1000);
  fMinTickDrift =
    clock_data.Time2Tick(clock_data.TriggerTime()); // this is the hardware trigger time
  fMaxTickDrift = fMinTickDrift + fDetectorWidthTicks + fEndTickPadding;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method. The goal is to produce a list of recob::Hit
/// objects which are a "clean" subset of all hits and which are believed to
/// be due to a neutrino interaction. It does this by considering input CosmicTag
/// objects, relating them to PFParticles/Tracks and removing the hits
/// associated to those objects which are believed to be Cosmic Rays.
///
void CRHitRemoval::produce(art::Event& evt)
{
  ++fNumEvent;

  // Start by looking up the original hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  evt.getByLabel(fHitProducerLabel, hitHandle);

  // If there are no hits then there should be no output
  if (!hitHandle.isValid()) return;

  // also get the associated wires
  // Note that we no longer assume RawDigit associations and will not preserve them if they exist
  // we assume they have been created by the same module as the hits
  art::FindOneP<recob::Wire> ChannelHitWires(hitHandle, evt, fHitProducerLabel);

  HitPtrVector ChHits;
  art::fill_ptr_vector(ChHits, hitHandle);

  // this object contains the hit collection
  // and its associations to wires and raw digits:
  recob::HitCollectionCreator hcol(evt, ChannelHitWires.isValid(), false);

  // Now recover thre remaining collections of objects in the event store that we need
  // Recover the PFParticles that are ultimately associated to tracks
  // This will be the key to recovering the hits
  art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
  evt.getByLabel(fPFParticleProducerLabel, pfParticleHandle);

  // Without a valid collection of PFParticles we can't do the hit removal
  if (!pfParticleHandle.isValid()) {
    copyAllHits(ChHits, ChannelHitWires, hcol);

    // put the hit collection and associations into the event
    hcol.put_into(evt);

    return;
  }

  // Recover the clusters so we can do associations to the hits
  // In theory the clusters come from the same producer as the PFParticles
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  evt.getByLabel(fPFParticleProducerLabel, clusterHandle);

  // If there are no clusters then something is really wrong
  if (!clusterHandle.isValid()) {
    copyAllHits(ChHits, ChannelHitWires, hcol);

    // put the hit collection and associations into the event
    hcol.put_into(evt);

    return;
  }

  // Build the list of cosmic tag producers, associations and thresholds
  std::vector<art::Handle<std::vector<anab::CosmicTag>>> cosmicHandleVec;
  std::vector<std::map<int, std::vector<art::Ptr<recob::Track>>>> cosmicToTrackVecMapVec;
  std::vector<std::map<int, std::vector<art::Ptr<recob::PFParticle>>>> trackToPFParticleVecMapVec;
  std::vector<double> thresholdVec;

  for (size_t idx = 0; idx != fCosmicProducerLabels.size(); idx++) {
    std::string& handleLabel(fCosmicProducerLabels[idx]);

    art::Handle<std::vector<anab::CosmicTag>> cosmicHandle;
    evt.getByLabel(handleLabel, cosmicHandle);

    if (cosmicHandle.isValid()) {
      // Look up the associations to tracks
      art::FindManyP<recob::Track> cosmicTrackAssns(cosmicHandle, evt, handleLabel);

      // Get the handle to the tracks
      art::Handle<std::vector<recob::Track>> trackHandle;
      evt.getByLabel(fTrackProducerLabels[idx], trackHandle);

      // This should be the case
      if (trackHandle.isValid()) {
        std::map<int, std::vector<art::Ptr<recob::Track>>> cosmicToTrackVecMap;
        std::map<int, std::vector<art::Ptr<recob::PFParticle>>> trackToPFParticleVecMap;

        // Loop through the cosmic tags
        for (size_t cosmicIdx = 0; cosmicIdx < cosmicHandle->size(); cosmicIdx++) {
          art::Ptr<anab::CosmicTag> cosmicTag(cosmicHandle, cosmicIdx);

          std::vector<art::Ptr<recob::Track>> cosmicToTrackVec =
            cosmicTrackAssns.at(cosmicTag.key());

          cosmicToTrackVecMap[cosmicTag.key()] = cosmicToTrackVec;

          art::FindManyP<recob::PFParticle> trackPFParticleAssns(
            trackHandle, evt, fAssnProducerLabels[idx]);

          for (auto& track : cosmicToTrackVec) {
            std::vector<art::Ptr<recob::PFParticle>> trackToPFParticleVec =
              trackPFParticleAssns.at(track.key());

            for (auto& pfParticle : trackToPFParticleVec)
              trackToPFParticleVecMap[track.key()].push_back(pfParticle);
          }
        }

        // Store these collections
        cosmicHandleVec.emplace_back(cosmicHandle);
        cosmicToTrackVecMapVec.emplace_back(cosmicToTrackVecMap);
        trackToPFParticleVecMapVec.emplace_back(trackToPFParticleVecMap);
        thresholdVec.push_back(fCosmicTagThresholds[idx]);
      }
    }
  }

  // No cosmic tags then nothing to do here
  if (cosmicHandleVec.empty()) {
    copyAllHits(ChHits, ChannelHitWires, hcol);

    // put the hit collection and associations into the event
    hcol.put_into(evt);

    return;
  }

  // Now we walk up the remaining list of associations needed to go from CR tags to hits
  // We'll need to go from tracks to PFParticles
  // From PFParticles we go to clusters
  art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleProducerLabel);

  // Likewise, recover the collection of associations to hits
  art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fPFParticleProducerLabel);

  // No point double counting hits
  std::set<const recob::PFParticle*> taggedSet;

  // Start the identification of hits to remove. The outer loop is over the various producers of
  // the CosmicTag objects we're examininig
  for (size_t idx = 0; idx != cosmicHandleVec.size(); idx++) {
    // Obviously, dereference the handle and associations
    const art::Handle<std::vector<anab::CosmicTag>>& cosmicHandle(cosmicHandleVec[idx]);
    const std::map<int, std::vector<art::Ptr<recob::Track>>>& crTagTrackVec(
      cosmicToTrackVecMapVec[idx]);
    const std::map<int, std::vector<art::Ptr<recob::PFParticle>>>& trackToPFParticleVecMap(
      trackToPFParticleVecMapVec[idx]);

    for (size_t crIdx = 0; crIdx != cosmicHandle->size(); crIdx++) {
      art::Ptr<anab::CosmicTag> cosmicTag(cosmicHandle, crIdx);

      // If this was tagged as a CR muon then we have work to do!
      if (cosmicTag->CosmicScore() > thresholdVec[idx]) {
        // Recover the associated track
        std::vector<art::Ptr<recob::Track>> trackVec(crTagTrackVec.at(cosmicTag.key()));

        // Loop over the tracks (almost always only 1)
        for (const auto& track : trackVec) {
          std::map<int, std::vector<art::Ptr<recob::PFParticle>>>::const_iterator
            trackToPFParticleVecItr = trackToPFParticleVecMap.find(track.key());

          if (trackToPFParticleVecItr != trackToPFParticleVecMap.end()) {
            // Loop through them
            for (const auto& pfParticlePtr : trackToPFParticleVecItr->second) {
              // Get bare pointer
              const recob::PFParticle* pfParticle = pfParticlePtr.get();

              // A cosmic ray must be a primary (by fiat)
              while (!pfParticle->IsPrimary())
                pfParticle =
                  art::Ptr<recob::PFParticle>(pfParticleHandle, pfParticle->Parent()).get();

              // Add to our list of tagged PFParticles
              taggedSet.insert(pfParticle);
            }
          }
        }
      }
    }
  }

  // If no PFParticles have been tagged then nothing to do
  if (!taggedSet.empty()) {
    // This may all seem backwards... but what you want to do is remove all the hits which are associated to tagged
    // cosmic ray tracks/PFParticles and what you want to leave is all other hits. This includes hits on other PFParticles
    // as well as unassociated (as yet) hits.
    // One also needs to deal with a slight complication with hits in 2D which are shared between untagged PFParticles and
    // tagged ones... we should leave these in.

    // All that is left is to go through the PFParticles which have not been tagged and collect up all their hits
    // to output to our new collection
    // Note that this SHOULD take care of the case of shared 2D hits automagically since if the PFParticle has not
    // been tagged and it shares hits we'll pick those up here.
    HitPtrVector taggedHits;
    HitPtrVector untaggedHits;

    // Loop through the PFParticles and build out the list of hits on untagged PFParticle trees
    for (const auto& pfParticle : *pfParticleHandle) {
      // Start with only primaries
      if (!pfParticle.IsPrimary()) continue;

      // Temporary container for these hits
      HitPtrVector tempHits;

      // Find the hits associated to this untagged PFParticle
      collectPFParticleHits(&pfParticle, pfParticleHandle, clusterAssns, clusterHitAssns, tempHits);

      // One more possible chance at identifying tagged hits...
      // Check these hits to see if any lie outside time window
      bool goodHits(true);

      if (taggedSet.find(&pfParticle) != taggedSet.end()) goodHits = false;

      if (goodHits) {
        int nOutOfTime(0);

        for (const auto& hit : tempHits) {
          // Check on out of time hits
          if (hit->PeakTimeMinusRMS() < fMinTickDrift || hit->PeakTimePlusRMS() > fMaxTickDrift)
            nOutOfTime++;

          if (nOutOfTime > fMaxOutOfTime) {
            goodHits = false;
            break;
          }
        }
      }

      if (goodHits)
        std::copy(tempHits.begin(), tempHits.end(), std::back_inserter(untaggedHits));
      else
        std::copy(tempHits.begin(), tempHits.end(), std::back_inserter(taggedHits));
    }

    // First task - remove hits from the tagged hit collection that are inthe untagged hits (shared hits)
    FilterHits(taggedHits, untaggedHits);

    // Now filter the tagged hits from the total hit collection
    FilterHits(ChHits, taggedHits);
  }

  // Copy our new hit collection to the output
  copyInTimeHits(ChHits, ChannelHitWires, hcol);

  // put the hit collection and associations into the event
  hcol.put_into(evt);
}

//----------------------------------------------------------------------------
/// Find all hits in PFParticle hierarchy
///
/// Arguments:
///
/// pfParticle - the top level PFParticle to have hits removed
/// pfParticleHandle - handle to the PFParticle objects
/// partToClusAssns - list of PFParticle to Cluster associations
/// clusToHitAssns - list of Cluster to Hit associations
/// hitVec - the current list of hits
///
/// This recursively called method will remove all hits associated to an input
/// PFParticle and, in addition, will call itself for all daughters of the input
/// PFParticle
///
void CRHitRemoval::collectPFParticleHits(
  const recob::PFParticle* pfParticle,
  const art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
  const art::FindManyP<recob::Cluster>& partToClusAssns,
  const art::FindManyP<recob::Hit>& clusToHitAssns,
  HitPtrVector& hitVec)
{
  // Recover the clusters associated to the input PFParticle
  std::vector<art::Ptr<recob::Cluster>> clusterVec = partToClusAssns.at(pfParticle->Self());

  int minTick(fMinTickDrift);
  int maxTick(fMaxTickDrift);

  // Loop over the clusters and grab the associated hits
  for (const auto& cluster : clusterVec) {
    std::vector<art::Ptr<recob::Hit>> clusHitVec = clusToHitAssns.at(cluster.key());
    hitVec.insert(hitVec.end(), clusHitVec.begin(), clusHitVec.end());

    int minHitTick = (*std::min_element(clusHitVec.begin(),
                                        clusHitVec.end(),
                                        [](const auto& hit, const auto& min) {
                                          return hit->PeakTimeMinusRMS() < min->PeakTimeMinusRMS();
                                        }))
                       ->PeakTimeMinusRMS();
    int maxHitTick = (*std::max_element(clusHitVec.begin(),
                                        clusHitVec.end(),
                                        [](const auto& hit, const auto& max) {
                                          return hit->PeakTimePlusRMS() > max->PeakTimePlusRMS();
                                        }))
                       ->PeakTimePlusRMS();

    if (minHitTick < minTick) minTick = minHitTick;
    if (maxHitTick < maxTick) maxTick = maxHitTick;
  }

  // Loop over the daughters of this particle and remove their hits as well
  for (const auto& daughterId : pfParticle->Daughters()) {
    art::Ptr<recob::PFParticle> daughter(pfParticleHandle, daughterId);

    collectPFParticleHits(
      daughter.get(), pfParticleHandle, partToClusAssns, clusToHitAssns, hitVec);
  }

  return;
}

void CRHitRemoval::copyAllHits(std::vector<art::Ptr<recob::Hit>>& inputHits,
                               art::FindOneP<recob::Wire>& wireAssns,
                               recob::HitCollectionCreator& newHitCollection)
{
  for (const auto& hitPtr : inputHits) {
    art::Ptr<recob::Wire> wire = wireAssns.at(hitPtr.key());

    // just copy it
    newHitCollection.emplace_back(*hitPtr, wire);
  }

  return;
}

void CRHitRemoval::copyInTimeHits(std::vector<art::Ptr<recob::Hit>>& inputHits,
                                  art::FindOneP<recob::Wire>& wireAssns,
                                  recob::HitCollectionCreator& newHitCollection)
{
  for (const auto& hitPtr : inputHits) {
    // Check on out of time hits
    if (hitPtr->PeakTimeMinusRMS() < fMinTickDrift || hitPtr->PeakTimePlusRMS() > fMaxTickDrift)
      continue;

    art::Ptr<recob::Wire> wire = wireAssns.at(hitPtr.key());

    // just copy it
    newHitCollection.emplace_back(*hitPtr, wire);
  }

  return;
}

//----------------------------------------------------------------------------
// Filter a collection of hits (set difference).
// This function is copied from Track3DKalmanHit_module.cc
//
// Arguments:
//
// hits      - Hit collection from which hits should be removed.
// used_hits - Hits to remove.
//
void CRHitRemoval::FilterHits(HitPtrVector& hits, HitPtrVector& used_hits)
{
  if (used_hits.size() > 0) {
    // Make sure both hit collections are sorted.
    std::stable_sort(hits.begin(), hits.end());
    std::stable_sort(used_hits.begin(), used_hits.end());

    // Do set difference operation.
    std::vector<art::Ptr<recob::Hit>>::iterator it = std::set_difference(
      hits.begin(), hits.end(), used_hits.begin(), used_hits.end(), hits.begin());

    // Truncate hit collection.
    hits.erase(it, hits.end());
  }
}

//----------------------------------------------------------------------------
/// End job method.
void CRHitRemoval::endJob()
{
  double aveCRPerEvent = fNumEvent > 0 ? double(fNumCRRejects) / double(fNumEvent) : 0.;

  mf::LogInfo("CRHitRemoval") << "CRHitRemoval statistics:\n"
                              << "  Number of events = " << fNumEvent << "\n"
                              << "  Number of Cosmic Rays found = " << fNumCRRejects << ", "
                              << aveCRPerEvent << " average/event";
}
