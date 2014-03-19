// Some kinda description here, maybe
//
// It does optical stuff.



#include "art/Framework/Core/EDProducer.h"
#include "AnalysisBase/FlashMatch.h"
#include "RecoAlg/SpacePointAlg.h"

// ROOT includes.
#include <Rtypes.h>
#ifndef FlashClusterMatch_h
#define FlashClusterMatch_h 1


namespace recob{
  class OpFlash;
  class Cluster;
  class Hit;
}


namespace opdet {
  

  class FlashClusterMatch : public art::EDProducer{
  public:
    
    FlashClusterMatch(const fhicl::ParameterSet&);
    virtual ~FlashClusterMatch();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    
    void beginJob();
    
    
  private:
    trkf::SpacePointAlg       *  fSptalg;

    std::string fClusterModuleLabel;
    std::string fFlashModuleLabel;
    int         fMinSptsForOverlap;
  };

}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  FlashClusterMatch_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(FlashClusterMatch)

}//end namespace opdet
////////////////////////////////////////////////////////////////////////


// LArSoft includes
#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/OpFlash.h"
#include "RecoBase/SpacePoint.h"

// FMWK includes
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/LArProperties.h"


// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>



namespace opdet {

  //-------------------------------------------------
  
  FlashClusterMatch::FlashClusterMatch(fhicl::ParameterSet const& pset)
  {
    produces< std::vector<anab::FlashMatch> >();
    produces< art::Assns<recob::Cluster, anab::FlashMatch> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void FlashClusterMatch::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");   
    fFlashModuleLabel   = pset.get<std::string>("FlashModuleLabel");
    fMinSptsForOverlap  = pset.get<int>("MinSptsForOverlap");
    fSptalg             = new trkf::SpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"));

  }


  //-------------------------------------------------

  void FlashClusterMatch::beginJob()
  {
  }



  //-------------------------------------------------

  FlashClusterMatch::~FlashClusterMatch()
  {
  }





  //-------------------------------------------------


  void FlashClusterMatch::produce(art::Event& evt)
  {

    int n=0;


  top:
    n++;

    // DO NOT REMOVE THIS LINE:
    mf::LogWarning("RecoBaseDefaultCtor") << "using default Hit ctor - should only ever"
      					  << " be done when getting hits out of an event"
      					  << " not when trying to produce new hits to store"
      					  << " in the event";

    std::cerr<< " Warning : you have disabled the RecoBaseDefaultCtor message."  ;
    std::cerr<< "  Should only ever be done when trying to avoid messages when getting hits out of an event, not when trying to produce new hits to store in the event."<< std::endl;
    
    ++n;
    if(n<10) goto top;



    
    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        Flashes.push_back(flash);
      }

    // Read in clusters from the event
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);
    std::vector<art::Ptr<recob::Cluster> >  Clusters;
    for(unsigned int i=0; i < clusterh->size(); ++i)
      {
	art::Ptr<recob::Cluster> cluster(clusterh,i);
	Clusters.push_back(cluster);
      }
    
    // Pull associated hits from event
    art::FindManyP<recob::Hit> hits(clusterh, evt, fClusterModuleLabel);


    std::vector<std::vector<int> > SortedByViews(3);

    //    std::vector<int> MaxWire(Clusters.size(), 0);
    //    std::vector<int> MinWire(Clusters.size(), 10000);

    std::vector<int> MaxTime(Clusters.size(), 0);
    std::vector<int> MinTime(Clusters.size(), 10000);


    // Sort clusters by view
    for(size_t iClus=0; iClus!=Clusters.size(); ++iClus)
      {
	SortedByViews[Clusters.at(iClus)->View()].push_back(iClus);
	for(size_t iHit=0; iHit!=hits.at(iClus).size(); ++iHit)
	  {
	    double Time = hits.at(iClus).at(iHit)->PeakTime();
	    if(Time > MaxTime[iClus]) MaxTime[iClus] = Time;
	    if(Time < MinTime[iClus]) MinTime[iClus] = Time;
	    
	    //  Equivalent info for wires, maybe we want it later.
	    //	    int Wire = hits.at(iClus).at(iHit)->WireID().Wire;	 
	    //	    if(Wire < MinWire[iClus]) MinWire[iClus] = Wire;
	    //	    if(Wire > MaxWire[iClus]) MaxWire[iClus] = Wire;

	  }
      }

    // Loop over sets of 3 clusters
    for(size_t nU=0; nU!=SortedByViews[0].size(); ++nU)
      for(size_t nV=0; nV!=SortedByViews[1].size(); ++nV)
	for(size_t nW=0; nW!=SortedByViews[2].size(); ++nW)
	  {
	    int indexU = SortedByViews[0][nU];
	    int indexV = SortedByViews[1][nV];
	    int indexW = SortedByViews[2][nW];

	    bool NoOverlap = false;
	    
	    // Skip over clusters with no time overlap
	    for(size_t v=0; v!=3; ++v)
	      {
		int v1 = (v+1)%3;
		int v2 = (v+2)%3;
		
		if(MinTime[v] > std::min(MaxTime[v1],MaxTime[v2]))
		  NoOverlap = true;

		if(MaxTime[v] < std::max(MinTime[v1],MinTime[v2]))
		  NoOverlap = true;
	      }
	    
	    if(NoOverlap) continue;

	    // Prepare flattened vector for space pointery
	    art::PtrVector<recob::Hit>  FlatHits;

	    FlatHits.insert(FlatHits.begin(), hits.at(indexU).begin(), hits.at(indexU).end());
	    FlatHits.insert(FlatHits.begin(), hits.at(indexV).begin(), hits.at(indexV).end());
	    FlatHits.insert(FlatHits.begin(), hits.at(indexW).begin(), hits.at(indexW).end());

	    // Make the spacepoints
	    std::vector<recob::SpacePoint> spts;
	    fSptalg->makeSpacePoints(FlatHits, spts);

	    if(int(spts.size()) < fMinSptsForOverlap) continue;
 	    
	  }

    std::unique_ptr< std::vector<anab::FlashMatch> > flash_matches ( new std::vector<anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::Cluster, anab::FlashMatch > > assn_cluster( new art::Assns<recob::Cluster, anab::FlashMatch>);





    evt.put(std::move(flash_matches));
    evt.put(std::move(assn_cluster));

    
  }


}

