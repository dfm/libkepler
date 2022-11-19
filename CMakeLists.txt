cmake_minimum_required(VERSION 3.11)
project(Kepler VERSION 0.1.0 LANGUAGES CXX)

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

# Testing using Catch
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.2.0
)
FetchContent_MakeAvailable(catch)

# SIMD support from xsimd
FetchContent_Declare(
  xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG        9.0.1
)
FetchContent_MakeAvailable(xsimd)



# include(ExternalProject)
# ExternalProject_Add(xsimd
#   GIT_REPOSITORY    https://github.com/xtensor-stack/xsimd.git
#   GIT_TAG           9.0.1
#   PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/xsimd-9.0.1
#   CONFIGURE_COMMAND ""
#   BUILD_COMMAND     ${CMAKE_COMMAND} -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR> <SOURCE_DIR>
#   INSTALL_COMMAND   make install
# )
# include_directories(${CMAKE_CURRENT_BINARY_DIR}/xsimd-9.0.1/include)

# # include(GNUInstallDirs)
# # set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS ON)


# # Get the git hash & print status
# execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE PROJECT_GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
# message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
# message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# # Default to Release build type
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
# endif()
# message(STATUS "Build-type: ${CMAKE_BUILD_TYPE}")

# # Export the list of compile-commands into compile_commands.json
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
# add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/benchmarks)

# # add_executable(demo demo.cc)
# # target_compile_options(demo PUBLIC -O3 -Wall -pedantic -Wextra -Werror)
