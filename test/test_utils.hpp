#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <xsimd/xsimd.hpp>

using namespace Catch::Matchers;
namespace xs = xsimd;

template <typename T>
struct default_abs {
  constexpr static T value = 1e-13;
};

template <>
struct default_abs<float> {
  constexpr static float value = 5e-5;
};

template <typename T>
struct default_rel {
  constexpr static T value = 5e-4;
};

template <>
struct default_rel<float> {
  constexpr static float value = 5e-2;
};

template <typename T>
struct tolerance {
  constexpr static typename T::value_type abs = default_abs<typename T::value_type>::value;
  constexpr static typename T::value_type rel = default_rel<typename T::value_type>::value;
};
