
add_subdirectory(OpHitFinder)

# OpFlashFinder_module.cc is not built - it seems specific to uboone

art_make( 
          LIB_LIBRARIES 
                           lardata_RecoBase
                           larana_OpticalDetector_OpHitFinder
                           larsim_MCSTReco
                           lardata_MCBase
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larana_OpticalDetector_OpDigiProperties_service
                           lardata_Utilities
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_UTILITIES}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larana_OpticalDetector_OpDigiProperties_service
                           larana_OpticalDetector_OpHitFinder
                           larana_OpticalDetector
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_MCCheater_BackTracker_service
                           larevt_Filters
                           lardata_RecoObjects
                           lardata_RecoBase
                           lardata_MCBase
                           larreco_RecoAlg
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation
                           larsim_TriggerAlgo_TriggerAlgoMicroBoone_service
                           lardata_Utilities
                           lardata_RawData
                           lardata_AnalysisAlg
                           lardata_AnalysisBase
                           SeedService_service # artextensions
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_UTILITIES}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${ROOT_BASIC_LIB_LIST}
          SERVICE_LIBRARIES
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_UTILITIES}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${CLHEP}
                           ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
install_fw( LIST toyWaveform.txt
            SUBDIRNAME OpticalDetector )