include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)

set(KEPLER_TESTS
  test_math
  test_reduction
  test_refiners
  test_solve
  test_starters)

foreach(name ${KEPLER_TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE Catch2::Catch2WithMain)

  if(WIN32)
    # TODO(dfm): Add a Windows-specific compiler flags
  else()
    target_compile_options(${name} PRIVATE -O3 -Wall -pedantic -Wextra -Werror)
    if(HAS_NO_DEPRECATED_COPY)
      target_compile_options(${name} PRIVATE -Wno-deprecated-copy)
    endif()
    if(HAS_MARCH_NATIVE)
      target_compile_options(${name} PRIVATE -march=native)
    endif()
  endif()

  add_test(${name} ${name})

endforeach()

add_executable(get_simd_size get_simd_size.cpp)
if(WIN32)
  # TODO(dfm): Add a Windows-specific compiler flags
else()
  target_compile_options(get_simd_size PRIVATE -O3 -Wall -pedantic -Wextra -Werror)
  if(HAS_NO_DEPRECATED_COPY)
    target_compile_options(get_simd_size PRIVATE -Wno-deprecated-copy)
  endif()
  if(HAS_MARCH_NATIVE)
    target_compile_options(get_simd_size PRIVATE -march=native)
  endif()
endif()

