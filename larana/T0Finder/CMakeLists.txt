add_subdirectory(AssociationsTools)

art_make(MODULE_LIBRARIES
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ROOT::Tree
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_AnalysisBase
           lardataobj_RecoBase
           larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
           nusimdata_SimulationBase)

install_headers()
install_fhicl()
install_source()
install_scripts()
