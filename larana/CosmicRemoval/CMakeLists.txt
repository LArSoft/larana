add_subdirectory(TrackContainment)

art_make(
  LIB_LIBRARIES
    ${ART_UTILITIES}
    fhiclcpp::fhiclcpp fhiclcpp::types
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    cetlib_except
    larcorealg_Geometry
    lardataobj_AnalysisBase
    lardataobj_RecoBase
    larsim_PhotonPropagation_PhotonVisibilityService_service
    nusimdata_SimulationBase
  MODULE_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    messagefacility::MF_MessageLogger
    ROOT::Core
    ROOT::Physics
    ROOT::Tree
    larana_CosmicRemoval
    larana_CosmicRemoval_TrackContainment
    lardataalg_DetectorInfo
    larcorealg_Geometry
    lardata_ArtDataHelper
    lardataobj_AnalysisBase
    lardataobj_RecoBase
    larreco_RecoAlg_Cluster3DAlgs
  )

install_headers()
install_fhicl()
install_source()
