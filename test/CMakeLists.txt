# Testing using Catch
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.2.0
)
FetchContent_MakeAvailable(Catch2)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-deprecated-copy HAS_NO_DEPRECATED_COPY)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)

set(KEPLER_TESTS
  test_reduction
  test_starters)

foreach(name ${KEPLER_TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE Catch2::Catch2WithMain)

  if(WIN32)
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
