cmake_minimum_required(VERSION 3.25)
project(AtlasMetFramework LANGUAGES CXX)

find_package(AnalysisBase 25.2.117 EXACT CONFIG REQUIRED)
if(NOT TARGET atlas_met_maker_correction)
    add_subdirectory(../maker_adapter maker_adapter)
endif()

add_library(AtlasMetFramework SHARED src/bridge.cpp src/execute_maker.cpp src/maker_xaod.cpp src/maker_trace.cpp)
target_compile_features(AtlasMetFramework PRIVATE cxx_std_23)
target_compile_definitions(AtlasMetFramework PRIVATE XAOD_STANDALONE XAOD_ANALYSIS)
target_link_libraries(AtlasMetFramework PRIVATE
    atlas_met_maker_correction
    AnalysisBase::ColumnarModeLib AnalysisBase::MetAnalysisAlgorithmsLib
    AnalysisBase::METUtilitiesLib AnalysisBase::xAODRootAccess
    AnalysisBase::xAODMissingET AnalysisBase::xAODJet AnalysisBase::xAODCore
    AnalysisBase::xAODMuon AnalysisBase::xAODTracking ${CMAKE_DL_LIBS}
)
if(UNIX AND NOT APPLE)
    target_link_options(AtlasMetFramework PRIVATE -Wl,-z,defs)
endif()

option(ATLAS_MET_FULL_BRIDGE "Use shared soft, Final and significance nodes" OFF)
if(ATLAS_MET_FULL_BRIDGE)
    if(NOT TARGET met_soft)
        add_subdirectory(../soft soft EXCLUDE_FROM_ALL)
    endif()
    target_compile_definitions(AtlasMetFramework PRIVATE ATLAS_MET_FULL_BRIDGE)
    target_sources(AtlasMetFramework PRIVATE src/soft_xaod.cpp src/execute_builder.cpp
        src/execute_significance.cpp)
    target_link_libraries(AtlasMetFramework PRIVATE atlas::met_soft AnalysisBase::xAODCaloEvent
        AnalysisBase::xAODEgamma AnalysisBase::xAODEventInfo AnalysisBase::xAODTau)
endif()
