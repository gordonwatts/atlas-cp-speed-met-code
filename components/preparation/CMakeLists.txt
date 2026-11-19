cmake_minimum_required(VERSION 3.25)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(MetPreparation LANGUAGES CXX)
endif()

option(MET_PREPARATION_STOCK_CHECK "Build the stock nominal preparation checkers" OFF)

if(NOT TARGET correctionlib)
    find_package(correctionlib CONFIG REQUIRED)
endif()
if(NOT TARGET met_significance_numeric)
    add_subdirectory(../significance "${CMAKE_CURRENT_BINARY_DIR}/significance" EXCLUDE_FROM_ALL)
endif()

add_library(met_preparation_nodes STATIC preparation.cxx tau.cxx muon.cxx egamma.cxx)
add_library(atlas::met_preparation ALIAS met_preparation_nodes)
target_link_libraries(met_preparation_nodes PUBLIC correctionlib PRIVATE met_significance_numeric)
target_include_directories(met_preparation_nodes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(met_preparation_nodes PUBLIC cxx_std_20)
set_target_properties(met_preparation_nodes PROPERTIES POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(met_preparation_nodes PRIVATE -ffp-contract=off)
endif()

add_library(met_preparation_provider MODULE module.cxx)
target_link_libraries(met_preparation_provider PRIVATE met_preparation_nodes)
set_target_properties(met_preparation_provider PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

if(MET_PREPARATION_STOCK_CHECK)
    if(NOT TARGET AnalysisBase::TauAnalysisToolsLib)
        find_package(AnalysisBase 25.2.117 EXACT CONFIG REQUIRED)
    endif()
    add_executable(met_tau_prepare_stock_check tau_stock.cxx)
    target_compile_features(met_tau_prepare_stock_check PRIVATE cxx_std_23)
    target_compile_definitions(met_tau_prepare_stock_check PRIVATE XAOD_STANDALONE XAOD_ANALYSIS)
    target_link_libraries(met_tau_prepare_stock_check PRIVATE correctionlib
        AnalysisBase::ColumnarModeLib AnalysisBase::xAODRootAccess AnalysisBase::TauAnalysisToolsLib)

    add_executable(met_muon_prepare_stock_check muon_stock.cxx)
    target_compile_features(met_muon_prepare_stock_check PRIVATE cxx_std_23)
    target_compile_definitions(met_muon_prepare_stock_check PRIVATE XAOD_STANDALONE XAOD_ANALYSIS)
    target_link_libraries(met_muon_prepare_stock_check PRIVATE correctionlib
        AnalysisBase::xAODRootAccess AnalysisBase::xAODMuon AnalysisBase::xAODEventInfo
        AnalysisBase::xAODTracking ROOT::TreePlayer)

    add_executable(met_egamma_prepare_stock_check egamma_stock.cxx egamma_boundaries.cxx)
    target_compile_features(met_egamma_prepare_stock_check PRIVATE cxx_std_23)
    target_compile_definitions(met_egamma_prepare_stock_check PRIVATE XAOD_STANDALONE XAOD_ANALYSIS)
    target_link_libraries(met_egamma_prepare_stock_check PRIVATE correctionlib
        AnalysisBase::xAODRootAccess AnalysisBase::xAODEgamma AnalysisBase::xAODCaloEvent
        AnalysisBase::xAODEventInfo AnalysisBase::xAODTracking AnalysisBase::egammaUtils
        AnalysisBase::IsolationCorrectionsLib ROOT::TreePlayer)
endif()
