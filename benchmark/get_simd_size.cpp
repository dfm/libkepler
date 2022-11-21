#include <iostream>
#include <xsimd/xsimd.hpp>

int main() {
  constexpr std::size_t float_size = xsimd::simd_type<float>::size;
  constexpr std::size_t double_size = xsimd::simd_type<double>::size;
  std::cout << "float size: " << float_size << std::endl;
  std::cout << "double size: " << double_size << std::endl;
  return 0;
}