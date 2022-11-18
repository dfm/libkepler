#ifndef KEPLER_MATH_CONSTANTS_HPP
#define KEPLER_MATH_CONSTANTS_HPP

#include <cstdint>

#include "./utils.hpp"

namespace kepler {
namespace constants {

// origin: xsimd/arch/xsimd_constants.hpp
/****************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 ****************************************************************************/

#define KB_DEFINE_CONSTANT(NAME, SINGLE, DOUBLE) \
  template <class T>                             \
  inline T NAME() noexcept {                     \
    return T(NAME<typename T::value_type>());    \
  }                                              \
  template <>                                    \
  inline float NAME<float>() noexcept {          \
    return bit_cast<float>((uint32_t)SINGLE);    \
  }                                              \
  template <>                                    \
  inline double NAME<double>() noexcept {        \
    return bit_cast<double>((uint64_t)DOUBLE);   \
  }

KB_DEFINE_CONSTANT(pi, 0x40490fdb, 0x400921fb54442d18)
KB_DEFINE_CONSTANT(twopi, 0x40c90fdb, 0x401921fb54442d18)
KB_DEFINE_CONSTANT(pio2, 0x3fc90fdb, 0x3ff921fb54442d18)
KB_DEFINE_CONSTANT(pio2_1, 0x3fc90f80, 0x3ff921fb54400000)
KB_DEFINE_CONSTANT(pio2_1t, 0x37354443, 0x3dd0b4611a626331)
KB_DEFINE_CONSTANT(pio2_2, 0x37354400, 0x3dd0b4611a600000)
KB_DEFINE_CONSTANT(pio2_2t, 0x2e85a308, 0x3ba3198a2e037073)
KB_DEFINE_CONSTANT(pio2_3, 0x2e85a300, 0x3ba3198a2e000000)
KB_DEFINE_CONSTANT(pio2_3t, 0x248d3132, 0x397b839a252049c1)
KB_DEFINE_CONSTANT(pio4, 0x3f490fdb, 0x3fe921fb54442d18)
KB_DEFINE_CONSTANT(twentypi, 0x427b53d1, 0x404f6a7a2955385e)
KB_DEFINE_CONSTANT(twoopi, 0x3f22f983, 0x3fe45f306dc9c883)
KB_DEFINE_CONSTANT(mediumpi, 0x43490fdb, 0x412921fb54442d18)

KB_DEFINE_CONSTANT(sixth, 0x3e2aaaab, 0x3fc5555555555555)
KB_DEFINE_CONSTANT(twentieth, 0x3d2aaaab, 0x3fa5555555555555)

KB_DEFINE_CONSTANT(nijenhuis_factor1, 0x40f4da39, 0x401e9b471164c596)
KB_DEFINE_CONSTANT(nijenhuis_factor2, 0x3fa6450f, 0x3ff4c8a1d518acbd)

#undef KB_DEFINE_CONSTANT

}  // namespace constants
}  // namespace kepler

#endif
