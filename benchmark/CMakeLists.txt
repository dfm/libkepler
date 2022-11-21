include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)

add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE Catch2::Catch2WithMain)

if(WIN32)
  # TODO(dfm): Add a Windows-specific compiler flags
else()
  target_compile_options(benchmark PRIVATE -O3 -Wall -pedantic -Wextra -Werror)
  if(HAS_NO_DEPRECATED_COPY)
    target_compile_options(benchmark PRIVATE -Wno-deprecated-copy)
  endif()
  if(HAS_MARCH_NATIVE)
    target_compile_options(benchmark PRIVATE -march=native)
  endif()
endif()
