# Set cmake version requirement
cmake_minimum_required(VERSION 3.9)

# Set the name of the project
project(xGEM2MT VERSION 1.0.1 LANGUAGES CXX)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MACOSX_RPATH ON)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Include this module to enable testing with CTest and dashboard submissions to CDash
include(CTest)

# Use ccache to speed up repeated compilations
include(CCache)

# Ensure proper configuration if in a conda environment
include(CondaAware)

# Set the version of the project
set(XGEM2MT_VERSION_MAJOR "1")
set(XGEM2MT_VERSION_MINOR "0")
set(XGEM2MT_VERSION_MICRO "1")
set(XGEM2MT_VERSION "${XGEM2MT_VERSION_MAJOR}.${XGEM2MT_VERSION_MINOR}.${XGEM2MT_VERSION_MICRO}")

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Define which components of Reaktoro to build
option(XGEM2MT_BUILD_ALL    "Build everything." OFF)
option(XGEM2MT_BUILD_PYTHON "Build the python bindings." ON)
option(XGEM2MT_BUILD_DEMOS  "Build demos." ON)
option(XGEM2MT_BUILD_DOCS   "Build documentation." OFF)
option(XGEM2MT_BUILD_TESTS  "Build tests." OFF)
option(XGEM2MT_USE_THERMOFUN "Use ThermoFun" ON)

# Modify the XGEM2MT_BUILD_X variables according to XGEM2MT_BUILD_ALL
if(XGEM2MT_BUILD_ALL MATCHES ON)
    set(XGEM2MT_BUILD_DEMOS  ON)
    set(XGEM2MT_BUILD_DOCS   ON)
    set(XGEM2MT_BUILD_PYTHON ON)
    set(XGEM2MT_BUILD_TESTS  ON)
endif()

# Set the CTest option BUILD_TESTING to the value of XGEM2MT_BUILD_TESTS
set(BUILD_TESTING XGEM2MT_BUILD_TESTS)

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_STATIC_LIBS "Build static libraries." ON)

# Set libraries to be compiled with position independent code mode (i.e., fPIC option in GNU compilers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Cannot currently produce a dynamic library using MSVC, only a static library
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(XGEM2MT_BUILD_SHARED_LIBS OFF)
endif()

# Used into conda only
if(DEFINED ENV{CONDA_PREFIX})
   option(USE_SPDLOG_PRECOMPILED "Use spdlog in compiled version" ON)
else()
   option(USE_SPDLOG_PRECOMPILED "Use spdlog in compiled version" OFF)
endif()

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# test thread safe
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=thread")

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
   add_compile_options( -pthread -Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   add_compile_options(-Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Intel compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    add_compile_options(-Wall -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /MP4
        /EHsc
        /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
        /DNOMINMAX
    )
endif()


# Set some necessary definitions
add_definitions(-DNODEARRAYLEVEL)
add_definitions(-DIPMGEMPLUGIN)

if(XGEM2MT_USE_THERMOFUN MATCHES ON)
  add_definitions(-DUSE_THERMOFUN)
  #add_definitions(-DUSE_THERMO_LOG)
endif()

# Set the project directory names
set(XGEM2MT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/xGEM2MT)
set(XGEM2MT_DEMOS_DIR  ${CMAKE_SOURCE_DIR}/demos)
set(XGEM2MT_DOCS_DIR   ${CMAKE_SOURCE_DIR}/docs)
set(XGEM2MT_TESTS_DIR  ${CMAKE_SOURCE_DIR}/tests)

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
#include_directories(${EIGEN3_INCLUDE_DIR})

# Find all dependencies
include(xGEM2MTFindDeps)

# Build the C++ library and the Python bindings
add_subdirectory(xGEM2MT)

# Build the demos
if(XGEM2MT_BUILD_DEMOS)
    add_subdirectory(demos)
endif()

# Build the docs
if(XGEM2MT_BUILD_DOCS)
    add_subdirectory(docs)
endif()

# Build the tests
if(XGEM2MT_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# Build python wraper
if(XGEM2MT_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# Install the cmake config files that permit users to use find_package(ChemicalFun)
include(xGEM2MTInstallCMakeConfigFiles)

