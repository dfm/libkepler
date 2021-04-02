set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)

file(GLOB_RECURSE all_sources RELATIVE ${CMAKE_CURRENT_LIST_DIR} *.cc)

add_custom_target(results_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/results)

foreach(source ${all_sources})
  get_filename_component(source_name ${source} NAME_WE)
  add_executable(${source_name} ${source})
  if(WIN32)
    target_compile_options(${source_name} PUBLIC /Wall /bigobj)
  else()
    target_compile_options(${source_name} PUBLIC -Wall -pedantic -Wextra -Werror -ffast-math)
    if(HAS_NO_DEPRECATED_COPY)
      target_compile_options(${source_name} PUBLIC -Wno-deprecated-copy)
    endif()
  endif()
#   set_property(TARGET ${source_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
#   add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})

  add_custom_target(
    ${source_name}_output ALL
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${source_name}
    DEPENDS ${source_name} results_directory
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/results/${source_name}.csv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results
  )

  message(STATUS "Adding benchmark: ${source_name}")
endforeach()
