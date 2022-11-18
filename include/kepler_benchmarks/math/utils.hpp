#ifndef KEPLER_MATH_UTILS_HPP
#define KEPLER_MATH_UTILS_HPP

#include <cstring>

namespace kepler {

template <class To, class From>
inline To bit_cast(From val) noexcept {
  static_assert(sizeof(From) == sizeof(To), "casting between compatible layout");
  To res;
  std::memcpy(&res, &val, sizeof(val));
  return res;
}

}  // namespace kepler

#endif
