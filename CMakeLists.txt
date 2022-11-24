cmake_minimum_required(VERSION 3.11)
project(Kepler VERSION 0.1.0 LANGUAGES CXX)

# Don't allow in-tree builds
if (CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
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

FetchContent_Declare(
  xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG        9.0.1
)
FetchContent_MakeAvailable(xsimd)
include_directories(${xsimd_SOURCE_DIR}/include)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.2.0
  )
  FetchContent_MakeAvailable(Catch2)

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()

  OPTION(BUILD_KEPLER_BENCHMARK "Build the benchmark script" OFF)
  if(BUILD_KEPLER_BENCHMARK)
    add_subdirectory(benchmark)
  endif()
endif()

add_subdirectory(src/kepler/metal)
