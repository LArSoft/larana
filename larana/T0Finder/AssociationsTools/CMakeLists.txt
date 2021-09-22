cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    ${ART_PERSISTENCY_PROVENANCE}
    messagefacility::MF_MessageLogger
    ROOT::Core
    canvas::canvas
    larcorealg_Geometry
    larsim_MCCheater_BackTrackerService_service
)

install_headers()
install_source()
install_fhicl()
