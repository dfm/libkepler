#ifndef KEPLER_KEPLER_HPP
#define KEPLER_KEPLER_HPP

#ifdef __cpp_if_constexpr
#define KEPLER_IF_CONSTEXPR if constexpr
#endif
#if !defined(KEPLER_IF_CONSTEXPR) && __cplusplus >= 201703L
#define KEPLER_IF_CONSTEXPR if constexpr
#endif
#if !defined(KEPLER_IF_CONSTEXPR)
#define KEPLER_IF_CONSTEXPR if
#endif

#include "./constants.hpp"
#include "./reduction.hpp"
#include "./refiners.hpp"
#include "./solve.hpp"
#include "./starters.hpp"

#endif