#ifndef _KEPLER_BENCHMARKS_CUDA_SUPPORT_H_
#define _KEPLER_BENCHMARKS_CUDA_SUPPORT_H_

#ifdef __CUDACC__
#define INLINE_OR_DEVICE __host__ __device__

template <class T>
INLINE_OR_DEVICE void swap(T& a, T& b) {
  T c(a);
  a = b;
  b = c;
}

#else
#define INLINE_OR_DEVICE inline

template <class T>
inline void swap(T& a, T& b) {
  std::swap(a, b);
}

template <class T>
inline void sincos(const T& x, T* sx, T* cx) {
  *sx = sin(x);
  *cx = cos(x);
}
#endif

#endif
