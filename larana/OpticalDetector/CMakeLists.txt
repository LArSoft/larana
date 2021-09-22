add_subdirectory(OpHitFinder)

# OpFlashFinder_module.cc is not built - it seems specific to uboone

art_make(
  LIB_LIBRARIES
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    art_root_io::tfile_support
    CLHEP::CLHEP
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::RIO
    ROOT::Tree
    canvas::canvas
    cetlib_except::cetlib_except
    larana_OpticalDetector_OpDigiProperties_service
    larana_OpticalDetector_OpHitFinder
    larcorealg_Geometry
    lardataobj_RecoBase
    larsim_MCSTReco
    larsim_PhotonPropagation_PhotonVisibilityService_service
    nusimdata_SimulationBase
  MODULE_LIBRARIES
    art::Framework_Services_Registry
    larana_OpticalDetector
    lardataobj_RawData
    larsim_MCCheater_ParticleInventoryService_service
    larsim_Simulation
    nurandom_RandomUtils_NuRandomService_service
  SERVICE_LIBRARIES
    CLHEP::CLHEP
    messagefacility::MF_MessageLogger
    ROOT::Core
    ROOT::Hist
    larcorealg_Geometry
  )

install_headers()
install_fhicl()
install_source()
install_fw(LIST toyWaveform.txt SUBDIRNAME OpticalDetector)
