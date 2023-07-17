#include <nanobind/nanobind.h>

#include <kepler/kepler.hpp>

int add(int a, int b) { return a + b; }

NB_MODULE(_libkepler, m) { m.def("add", &add); }
