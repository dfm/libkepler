include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)

set(KEPLER_TESTS
  test_reduction
  test_refiners
  test_solve
  test_starters)

foreach(name ${KEPLER_TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE Catch2::Catch2WithMain)

  if(MSVC)
    target_compile_options(${name} PRIVATE /arch:AVX2)
  else()
    target_compile_options(${name} PRIVATE -O3 -Wall -pedantic -Wextra -Werror)
    if(HAS_NO_DEPRECATED_COPY)
      target_compile_options(${name} PRIVATE -Wno-deprecated-copy)
    endif()
    if(HAS_MARCH_NATIVE AND NOT CMAKE_CXX_FLAGS MATCHES "-march" AND NOT CMAKE_CXX_FLAGS MATCHES "-arch" AND NOT CMAKE_OSX_ARCHITECTURES)
      target_compile_options(${name} PRIVATE -march=native -mtune=native)
    endif()
  endif()

  add_test(${name} ${name})
endforeach()
