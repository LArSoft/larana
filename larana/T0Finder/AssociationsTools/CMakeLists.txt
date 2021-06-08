cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_PERSISTENCY_PROVENANCE}
    ${MF_MESSAGELOGGER}
    ROOT::Core
    canvas
    larcorealg_Geometry
    larsim_MCCheater_BackTrackerService_service
)

install_headers()
install_source()
install_fhicl()
