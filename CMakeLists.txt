# ======================================================================
#  larana main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larana/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.12.00 REQUIRED)
project(larana VERSION 10.01.06 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(art REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)

find_package(nug4 REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)

find_package(CLHEP COMPONENTS Random REQUIRED EXPORT)
find_package(Eigen3 REQUIRED)
find_package(ROOT COMPONENTS Core GenVector Hist MathCore Physics RIO TMVA Tree REQUIRED EXPORT)

find_package(larcore REQUIRED EXPORT)
find_package(larcorealg REQUIRED EXPORT)
find_package(larcoreobj REQUIRED EXPORT)
find_package(lardata REQUIRED EXPORT)
find_package(lardataalg REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(larreco REQUIRED EXPORT)
find_package(larsim REQUIRED EXPORT)

cet_cmake_module_directories(Modules BINARY)

# macros for dictionary and simple_plugin
include(BuildDictionary)
include(CetMake)
include(BuildPlugins)

# source
add_subdirectory(larana)

# tests
add_subdirectory(test)

# ups - table and config files

# packaging utility

cet_cmake_config()
