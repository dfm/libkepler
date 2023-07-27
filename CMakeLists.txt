cmake_minimum_required(VERSION 3.15)
project(Kepler LANGUAGES CXX)

# Don't allow in-tree builds
if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(
    FATAL_ERROR
    "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

set(KEPLER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versioning
file(STRINGS "${KEPLER_INCLUDE_DIR}/kepler/config.hpp" kepler_version_defines
  REGEX "#define KEPLER_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${kepler_version_defines})
  if(ver MATCHES "#define KEPLER_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
    set(KEPLER_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
  endif()
endforeach()

set(${PROJECT_NAME}_VERSION
  ${KEPLER_VERSION_MAJOR}.${KEPLER_VERSION_MINOR}.${KEPLER_VERSION_PATCH})
message(STATUS "libkepler v${${PROJECT_NAME}_VERSION}")

# We need modern C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Include the source for the libkepler project
include_directories(${CMAKE_CURRENT_LIST_DIR}/extern/vendor/xsimd-10.0.0)
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL "Kepler" AND NOT SKBUILD)
  include(CTest)
  option(KEPLER_BUILD_BENCHMARK "Build the benchmark script" OFF)

  if(BUILD_TESTING OR KEPLER_BUILD_BENCHMARK)
    include(FetchContent)

    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()

  if(KEPLER_BUILD_BENCHMARK)
    add_subdirectory(benchmark)
  endif()
endif()

# Build the Python extension -- this is really just for development purposes;
# specifically to provide the appropriate VSCode support via CMake. To actually
# build the Python extension for use or testing, it's better to use pip or build
# directly.
option(KEPLER_BUILD_PYTHON "Build the Python extension" OFF)

if(KEPLER_BUILD_PYTHON)
  include(FetchContent)

  FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v1.4.0
  )
  FetchContent_MakeAvailable(nanobind)

  # This is a hack that is needed since we have a check that the Python bindings
  # only get built by scikit-build-core by default.
  set(SKBUILD ON)

  add_subdirectory(python)
endif()
