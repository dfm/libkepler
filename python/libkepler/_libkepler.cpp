#include <nanobind/nanobind.h>

#include <iostream>
#include <kepler/kepler.hpp>
#include <string_view>
#include <xsimd/xsimd.hpp>

namespace specs {
// Ref:
// https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/64490578#64490578
template <typename T>
constexpr std::string_view type_name();

template <>
constexpr std::string_view type_name<void>() {
  return "void";
}

namespace detail {

using type_name_prober = void;

template <typename T>
constexpr std::string_view wrapped_type_name() {
#ifdef __clang__
  return __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
  return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
  return __FUNCSIG__;
#else
#error "Unsupported compiler"
#endif
}

constexpr std::size_t wrapped_type_name_prefix_length() {
  return wrapped_type_name<type_name_prober>().find(type_name<type_name_prober>());
}

constexpr std::size_t wrapped_type_name_suffix_length() {
  return wrapped_type_name<type_name_prober>().length() - wrapped_type_name_prefix_length() -
         type_name<type_name_prober>().length();
}

}  // namespace detail

template <typename T>
constexpr std::string_view type_name() {
  constexpr auto wrapped_name = detail::wrapped_type_name<T>();
  constexpr auto prefix_length = detail::wrapped_type_name_prefix_length();
  constexpr auto suffix_length = detail::wrapped_type_name_suffix_length();
  constexpr auto type_name_length = wrapped_name.length() - prefix_length - suffix_length;
  return wrapped_name.substr(prefix_length, type_name_length);
}

template <typename... Args>
void print_types(const xsimd::arch_list<Args...>&) {
  ((std::cout << type_name<Args>() << " "), ...);
  std::cout << std::endl;
}
}  // namespace specs

int print_system_specs() {
  std::cout << "Supported architectures: ";
  xsimd::supported_architectures archs;
  specs::print_types(archs);

  std::cout << "Best architecture: "
            << specs::type_name<xsimd::detail::best<xsimd::supported_architectures>::type>()
            << std::endl;

  constexpr std::size_t float_size = xsimd::simd_type<float>::size;
  constexpr std::size_t double_size = xsimd::simd_type<double>::size;
  std::cout << "SIMD size (float): " << float_size << std::endl;
  std::cout << "SIMD size (double): " << double_size << std::endl;
  return 0;
}

int add(int a, int b) { return a + b; }

NB_MODULE(_libkepler, m) {
  m.def("print_system_specs", &print_system_specs);
  m.def("add", &add);
}
