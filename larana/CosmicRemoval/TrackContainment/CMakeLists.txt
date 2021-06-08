art_make(
  LIB_LIBRARIES
    ${FHICLCPP}
    ROOT::Core
    cetlib_except
    larcorealg_Geometry
    lardataobj_AnalysisBase
    lardataobj_RecoBase
  )



install_headers()
install_source()
