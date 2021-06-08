# OpFlashFinder_module.cc is not built - it seems specific to uboone
art_make(
  LIB_LIBRARIES
    lardataobj_RecoBase
    larcorealg_Geometry
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib_except
    ROOT::Core
    ROOT::Hist
  )

install_headers()
install_source()
