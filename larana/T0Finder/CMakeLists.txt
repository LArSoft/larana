add_subdirectory(AssociationsTools)

art_make(MODULE_LIBRARIES
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${MF_MESSAGELOGGER}
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
