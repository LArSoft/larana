cet_make_library(SOURCE
  PhysdEdx.cxx
  Chi2PIDAlg.cxx
	MVAAlg.cxx
	PIDAAlg.cxx
  LikelihoodPIDAlg.cxx
  LIBRARIES
  PUBLIC
  larreco::Calorimetry
  larreco::RecoAlg
  lardataobj::AnalysisBase
  larcoreobj::SimpleTypesAndConstants
  canvas::canvas
  ROOT::GenVector
  ROOT::Hist
  ROOT::Physics
  ROOT::TMVA
  ROOT::MathMore
  PRIVATE
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::MathCore
  ROOT::RIO
  ROOT::Tree
)

cet_build_plugin(Chi2ParticleID art::EDProducer
  LIBRARIES PRIVATE
  larana::ParticleIdentification
  lardata::AssociationUtil
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MVAPID art::EDProducer
  LIBRARIES PRIVATE
  larana::ParticleIdentification
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Tree
)

cet_build_plugin(PIDAAnalyzer art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::ParticleIdentification
  larana::OpticalDetector
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Hist
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
