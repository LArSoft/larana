/**
 * \file PMTAna_module.cc
 *
 * \ingroup PMTAna
 *
 * \brief Class definition file of PMTAna
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PMTAna

@{*/

// ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/fwd.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft
#include "lardataobj/RawData/OpDetWaveform.h"

// STL
#include <functional>
#include <numeric>
#include <string>

// ROOT
#include <TTree.h>

// My modules
#include "OpHitFinder/AlgoFixedWindow.h"
#include "OpHitFinder/AlgoThreshold.h"
#include "OpHitFinder/PedAlgoEdges.h"
#include "OpHitFinder/PulseRecoManager.h"

namespace pmtana {

  /**
     \class PMTAna
     PMTAna module to copy LArSoft data contents into LArLight data formatted file
  */
  class PMTAna : public art::EDAnalyzer {

  public:
    /// Constructor
    PMTAna(const fhicl::ParameterSet&);

    /// Function to be called per event
    void analyze(const art::Event&);

  private:
    std::string _fifo_mod_name; ///< Input FIFOChannel producer name
    TTree* _tree;               ///< output data holder TTree

    PulseRecoManager _preco_man;
    AlgoThreshold _th_algo;
    AlgoFixedWindow _fw_algo;
    PedAlgoEdges _ped_algo;
  };

}

namespace pmtana {
  DEFINE_ART_MODULE(PMTAna)
}

namespace pmtana {

  //#######################################################################################################
  PMTAna::PMTAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), _preco_man(), _th_algo(), _fw_algo(), _ped_algo()
  //#######################################################################################################
  {

    // Obtain module names for input data
    _fifo_mod_name = pset.get<std::string>("fModName_FIFOChannel");

    // Next we make storage data class objects for those data types specified in fcl files.
    art::ServiceHandle<art::TFileService const> fileService;

    // Create TTree
    _tree = fileService->make<TTree>("pmt_tree", "Analysis Tree");

    //
    // Demonstration purpose ...
    //
    _preco_man.AddRecoAlgo(&_th_algo);
    _preco_man.AddRecoAlgo(&_fw_algo);
    _preco_man.SetDefaultPedAlgo(&_ped_algo);
  }

  //#######################################################################################################
  void PMTAna::analyze(const art::Event& evt)
  //#######################################################################################################
  {
    auto const& pmts = evt.getHandle<std::vector<raw::OpDetWaveform>>(_fifo_mod_name);
    if (!pmts) { return; }

    for (raw::OpDetWaveform const& waveform : *pmts) {
      _preco_man.Reconstruct(waveform);
    }
  }

}

/** @}*/ // end of PMTAna group
