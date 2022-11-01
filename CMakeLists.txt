cmake_minimum_required(VERSION 3.5)
project(KeplerBenchmarks VERSION 0.1.0 LANGUAGES CXX)

if (CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

include(ExternalProject)
ExternalProject_Add(xsimd
  GIT_REPOSITORY    https://github.com/xtensor-stack/xsimd.git
  GIT_TAG           9.0.1
  PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/xsimd-9.0.1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ${CMAKE_COMMAND} -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR> <SOURCE_DIR> 
  INSTALL_COMMAND   make install
)

# include(GNUInstallDirs)
# set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Get the git hash & print status
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE PROJECT_GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message(STATUS "Build-type: ${CMAKE_BUILD_TYPE}")

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/benchmarks)
