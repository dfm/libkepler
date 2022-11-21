include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)

add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE Catch2::Catch2WithMain)

if(MSVC)
  # TODO(dfm): Add a Windows-specific compiler flags
  target_compile_options(benchmark PRIVATE /arch:AVX2 /O2)
else()
  target_compile_options(benchmark PRIVATE -O3 -Wall -pedantic -Wextra -Werror)
  if(HAS_NO_DEPRECATED_COPY)
    target_compile_options(benchmark PRIVATE -Wno-deprecated-copy)
  endif()
  if(HAS_MARCH_NATIVE AND NOT CMAKE_CXX_FLAGS MATCHES "-march" AND NOT CMAKE_CXX_FLAGS MATCHES "-arch" AND NOT CMAKE_OSX_ARCHITECTURES)
    target_compile_options(benchmark PRIVATE -march=native -mtune=native)
  endif()
endif()
