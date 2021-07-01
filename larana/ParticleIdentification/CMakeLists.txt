art_make(
  LIB_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Provenance
    art_root_io::TFileService_service
    art_root_io::tfile_support
    art_Utilities
    fhiclcpp::fhiclcpp fhiclcpp::types
    messagefacility::MF_MessageLogger
    ${CLHEP_LIBRARIES}
    ROOT::Core
    ROOT::Hist
    ROOT::MathCore
    ROOT::Physics
    ROOT::RIO
    ROOT::TMVA
    ROOT::Tree
    canvas
    cetlib
    cetlib_except
    larcorealg_Geometry
    lardataobj_AnalysisBase
    larreco_Calorimetry
    nusimdata_SimulationBase
  MODULE_LIBRARIES
   larana_ParticleIdentification
  )

install_headers()
install_fhicl()

add_subdirectory(scripts)
