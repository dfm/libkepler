cmake_minimum_required(VERSION 3.11)
project(Kepler VERSION 0.1.0 LANGUAGES CXX)

# Don't allow in-tree builds
if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(
    FATAL_ERROR
    "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

# We need modern C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Add external dependencies
include(FetchContent)

# https://cliutils.gitlab.io/modern-cmake/chapters/projects/fetch.html
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  macro(FetchContent_MakeAvailable NAME)
    FetchContent_GetProperties(${NAME})

    if(NOT ${NAME}_POPULATED)
      FetchContent_Populate(${NAME})
      add_subdirectory(${${NAME}_SOURCE_DIR} ${${NAME}_BINARY_DIR})
    endif()
  endmacro()
endif()

# xsimd -- Cross-platform SIMD support
FetchContent_Declare(
  xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG 10.0.0
)

# Instead of
# > FetchContent_MakeAvailable(xsimd)
# we use the following, so that the xsimd headers don't get installed
# ref: https://stackoverflow.com/questions/65527126/disable-install-for-fetchcontent
FetchContent_GetProperties(xsimd)

if(NOT xsimd_POPULATED)
  FetchContent_Populate(xsimd)
  add_subdirectory(${xsimd_SOURCE_DIR} ${xsimd_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

include_directories(${xsimd_SOURCE_DIR}/include)

# Include the source for the libkepler project
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL "Kepler" AND NOT SKBUILD)
  include(CTest)
  option(KEPLER_BUILD_BENCHMARK "Build the benchmark script" OFF)

  if(BUILD_TESTING OR KEPLER_BUILD_BENCHMARK)
    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG v3.2.0
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