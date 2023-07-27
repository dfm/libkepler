cmake_minimum_required(VERSION 3.15...3.26)

if(NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# We need modern C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Include the source for the libkepler project
include_directories(${CMAKE_CURRENT_LIST_DIR}/../extern/vendor/xsimd-10.0.0)
include_directories(${CMAKE_CURRENT_LIST_DIR}/../include)

# Below is all just the standard scikit-build-core boilerplate.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(
  _libkepler
  STABLE_ABI
  NB_STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/libkepler/_libkepler.cpp
)
install(TARGETS _libkepler LIBRARY DESTINATION libkepler)
