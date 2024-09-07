#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <xsimd/xsimd.hpp>

#include "kepler/kepler/starters.hpp"

using namespace Catch::Matchers;
namespace xs = xsimd;

template <typename T>
struct default_abs {
  constexpr static T value = 1e-13;
};

template <>
struct default_abs<float> {
  constexpr static float value = 5e-5f;
};

template <typename T>
struct default_rel {
  constexpr static T value = 5e-4;
};

template <>
struct default_rel<float> {
  constexpr static float value = 5e-2f;
};

template <typename T>
struct tolerance {
  constexpr static typename T::value_type abs = default_abs<typename T::value_type>::value;
  constexpr static typename T::value_type rel = default_rel<typename T::value_type>::value;
};

template <typename R, typename S = kepler::starters::basic<typename R::value_type>>
struct SolveTestCase {
  typedef typename R::value_type value_type;
  typedef R refiner_type;
  typedef S starter_type;
};
