cmake_minimum_required(VERSION 3.23...3.30)
project(Kepler VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(
    FATAL_ERROR
    "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Set the default build type to Release" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

include(FetchContent)
FetchContent_Declare(
  xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG 13.0.0
)
FetchContent_MakeAvailable(xsimd)

add_library(kepler SHARED src/kepler.cpp)
target_include_directories(
  kepler PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_include_directories(kepler PRIVATE ${xsimd_SOURCE_DIR}/include)

include(GNUInstallDirs)
install(TARGETS kepler PUBLIC_HEADER)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
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
