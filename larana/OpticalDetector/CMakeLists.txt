add_subdirectory(OpHitFinder)

# OpFlashFinder_module.cc is not built - it seems specific to uboone

art_make(
  LIB_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${CLHEP_LIBRARIES}
    ${FHICLCPP}
    ${MF_MESSAGELOGGER}
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::RIO
    ROOT::Tree
    canvas
    cetlib_except
    larana_OpticalDetector_OpDigiProperties_service
    larana_OpticalDetector_OpHitFinder
    larcorealg_Geometry
    lardataobj_RecoBase
    larsim_MCSTReco
    larsim_PhotonPropagation_PhotonVisibilityService_service
    nusimdata_SimulationBase
  MODULE_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    larana_OpticalDetector
    lardataobj_RawData
    larsim_MCCheater_ParticleInventoryService_service
    larsim_Simulation
    nurandom_RandomUtils_NuRandomService_service
  SERVICE_LIBRARIES
    ${CLHEP_LIBRARIES}
    ${MF_MESSAGELOGGER}
    ROOT::Core
    ROOT::Hist
    larcorealg_Geometry
  )

install_headers()
install_fhicl()
install_fw(LIST toyWaveform.txt SUBDIRNAME OpticalDetector)
