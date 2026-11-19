cmake_minimum_required(VERSION 3.25)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(MetSignificance LANGUAGES CXX)
endif()

option(MET_SIGNIFICANCE_STOCK_CHECK "Build the release comparison with AnalysisBase" OFF)

add_library(met_significance_numeric STATIC resolution.cxx payload.cxx reduction.cxx)
target_include_directories(met_significance_numeric PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(met_significance_numeric PUBLIC cxx_std_20)
set_target_properties(met_significance_numeric PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(met_significance_numeric PRIVATE -ffp-contract=off)
endif()

if(NOT TARGET correctionlib)
    find_package(correctionlib CONFIG REQUIRED)
endif()
add_library(met_significance_nodes STATIC plugin.cxx)
target_link_libraries(met_significance_nodes PUBLIC met_significance_numeric correctionlib)
set_target_properties(met_significance_nodes PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(met_significance_provider MODULE module.cxx)
target_link_libraries(met_significance_provider PRIVATE met_significance_nodes)

add_executable(met_significance_smoke provider_smoke.cxx)
target_link_libraries(met_significance_smoke PRIVATE met_significance_numeric correctionlib)

if(MET_SIGNIFICANCE_STOCK_CHECK)
    if(NOT TARGET AnalysisBase::METUtilitiesLib)
        find_package(AnalysisBase 25.2.117 EXACT CONFIG REQUIRED)
    endif()
    add_executable(met_significance_stock_compare stock_compare.cxx)
    target_compile_features(met_significance_stock_compare PRIVATE cxx_std_23)
    target_compile_definitions(met_significance_stock_compare PRIVATE XAOD_STANDALONE XAOD_ANALYSIS)
    target_link_libraries(met_significance_stock_compare PRIVATE
        met_significance_numeric correctionlib
        AnalysisBase::ColumnarModeLib AnalysisBase::xAODRootAccess AnalysisBase::METUtilitiesLib
        AnalysisBase::JetCalibToolsLib AnalysisBase::ElectronPhotonFourMomentumCorrectionLib
        AnalysisBase::MuonMomentumCorrectionsLib AnalysisBase::TauAnalysisToolsLib
    )
endif()
