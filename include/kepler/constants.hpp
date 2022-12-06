#ifndef KEPLER_CONSTANTS_HPP
#define KEPLER_CONSTANTS_HPP

#include <cstdint>

#include "utils.hpp"

namespace kepler {
namespace constants {

#define KEPLER_DEFINE_CONSTANT(NAME, SINGLE, DOUBLE) \
  template <typename T>                              \
  inline T NAME() noexcept {                         \
    return T(NAME<typename T::value_type>());        \
  }                                                  \
  template <>                                        \
  inline float NAME<float>() noexcept {              \
    return bit_cast<float>((uint32_t)SINGLE);        \
  }                                                  \
  template <>                                        \
  inline double NAME<double>() noexcept {            \
    return bit_cast<double>((uint64_t)DOUBLE);       \
  }

KEPLER_DEFINE_CONSTANT(pi, 0x40490fdb, 0x400921fb54442d18)
KEPLER_DEFINE_CONSTANT(twopi, 0x40c90fdb, 0x401921fb54442d18)
KEPLER_DEFINE_CONSTANT(pio2, 0x3fc90fdb, 0x3ff921fb54442d18)
KEPLER_DEFINE_CONSTANT(pio3, 0x3f860a92, 0x3ff0c152382d7365)
KEPLER_DEFINE_CONSTANT(pio4, 0x3f490fdb, 0x3fe921fb54442d18)
KEPLER_DEFINE_CONSTANT(pio6, 0x3f060a92, 0x3fe0c152382d7365)
KEPLER_DEFINE_CONSTANT(pio12, 0x3e860a92, 0x3fd0c152382d7365)

KEPLER_DEFINE_CONSTANT(twopio3, 0x40060a92, 0x4000c152382d7365)
KEPLER_DEFINE_CONSTANT(threepio4, 0x4016cbe4, 0x4002d97c7f3321d2)

KEPLER_DEFINE_CONSTANT(fivepio6, 0x40278d36, 0x4004f1a6c638d03f)
KEPLER_DEFINE_CONSTANT(fivepio12, 0x3fa78d36, 0x3ff4f1a6c638d03f)
KEPLER_DEFINE_CONSTANT(sevenpio12, 0x3fea927f, 0x3ffd524fe24f89f1)
KEPLER_DEFINE_CONSTANT(elevenpio12, 0x40384e88, 0x400709d10d3e7eab)

// Limits for range reduction
KEPLER_DEFINE_CONSTANT(twentypi, 0x427b53d1, 0x404f6a7a2955385e)
KEPLER_DEFINE_CONSTANT(twoopi, 0x3f22f983, 0x3fe45f306dc9c883)
KEPLER_DEFINE_CONSTANT(mediumpi, 0x43490fdb, 0x412921fb54442d18)

KEPLER_DEFINE_CONSTANT(range_max, 0x4b800000, 0x4340000000000000)

// Higher precision digits of pi/2
KEPLER_DEFINE_CONSTANT(pio2_1, 0x3fc90f80, 0x3ff921fb54400000)
KEPLER_DEFINE_CONSTANT(pio2_1t, 0x37354443, 0x3dd0b4611a626331)
KEPLER_DEFINE_CONSTANT(pio2_2, 0x37354400, 0x3dd0b4611a600000)
KEPLER_DEFINE_CONSTANT(pio2_2t, 0x2e85a308, 0x3ba3198a2e037073)
KEPLER_DEFINE_CONSTANT(pio2_3, 0x2e85a300, 0x3ba3198a2e000000)
KEPLER_DEFINE_CONSTANT(pio2_3t, 0x248d3132, 0x397b839a252049c1)

KEPLER_DEFINE_CONSTANT(sixth, 0x3e2aaaab, 0x3fc5555555555555)
KEPLER_DEFINE_CONSTANT(twentieth, 0x3d4ccccd, 0x3fa999999999999a)

// Algorithm-specific constants
KEPLER_DEFINE_CONSTANT(markley_factor1, 0x40f4da39, 0x401e9b471164c596)
KEPLER_DEFINE_CONSTANT(markley_factor2, 0x3fa6450f, 0x3ff4c8a1d518acbd)

KEPLER_DEFINE_CONSTANT(rppb_g2s, 0x3e8483ee, 0x3fd0907dc1930690)
KEPLER_DEFINE_CONSTANT(rppb_g3s, 0x3f000000, 0x3fe0000000000000)
KEPLER_DEFINE_CONSTANT(rppb_g4s, 0x3f3504f3, 0x3fe6a09e667f3bcc)
KEPLER_DEFINE_CONSTANT(rppb_g5s, 0x3f5db3d7, 0x3febb67ae8584caa)
KEPLER_DEFINE_CONSTANT(rppb_g6s, 0x3f7746ea, 0x3feee8dd4748bf15)

// Coefficients for the series expansion of sin(x)
KEPLER_DEFINE_CONSTANT(shortsin1, 0x3e2aaaab, 0x3fc5555555555555)
KEPLER_DEFINE_CONSTANT(shortsin2, 0x3c088889, 0x3f81111111111111)
KEPLER_DEFINE_CONSTANT(shortsin3, 0x39500d01, 0x3f2a01a01a01a01a)
KEPLER_DEFINE_CONSTANT(shortsin4, 0x3638ef1d, 0x3ec71de3a556c734)
KEPLER_DEFINE_CONSTANT(shortsin5, 0x32d7322b, 0x3e5ae64567f544e4)
KEPLER_DEFINE_CONSTANT(shortsin6, 0x2f309231, 0x3de6124613a86d09)
KEPLER_DEFINE_CONSTANT(shortsin7, 0x2b573f9f, 0x3d6ae7f3e733b81f)

// Factorials for Householder's method
KEPLER_DEFINE_CONSTANT(hh2, 0x3f000000, 0x3fe0000000000000)
KEPLER_DEFINE_CONSTANT(hh3, 0x3e2aaaab, 0x3fc5555555555555)
KEPLER_DEFINE_CONSTANT(hh4, 0x3d2aaaab, 0x3fa5555555555555)
KEPLER_DEFINE_CONSTANT(hh5, 0x3c088889, 0x3f81111111111111)
KEPLER_DEFINE_CONSTANT(hh6, 0x3ab60b61, 0x3f56c16c16c16c17)
KEPLER_DEFINE_CONSTANT(hh7, 0x39500d01, 0x3f2a01a01a01a01a)

#undef KEPLER_DEFINE_CONSTANT

/* Note to self: to generate HEX constants in Python:

>>> import struct
>>> print(hex(struct.unpack('!L', struct.pack('!f', v))[0]))
>>> print(hex(struct.unpack('!Q', struct.pack('!d', v))[0]))
*/

}  // namespace constants
}  // namespace kepler

#endif
