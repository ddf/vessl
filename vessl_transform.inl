#pragma once

namespace vessl
{
template <typename T>
VESSL_INLINE transform::complex<T>& transform::complex<T>::scale(T scalar)
{
  r *= scalar;
  i *= scalar;
  return *this;
}

template <typename T>
VESSL_INLINE transform::complex<T>& transform::complex<T>::add(const complex &other)
{
  r += other.r;
  i += other.i;
  return *this;
}

template <typename T>
VESSL_INLINE transform::complex<T>& transform::complex<T>::subtract(const complex &other)
{
  r -= other.r;
  i -= other.i;
  return *this;
}

template <typename T>
VESSL_INLINE transform::complex<T>& transform::complex<T>::multiply(const complex &other)
{
  T rr = r*other.r - i*other.i;
  T ii = r*other.i + i*other.r;
  r = rr;
  i = ii;
  return *this;
}

template <typename T>
VESSL_INLINE transform::complex<T>& transform::complex<T>::operator=(const complex &other)
{
  r = other.r;
  i = other.i;
  return *this;
}

template <typename T>
VESSL_INLINE void transform::complex<T>::set_complex(T real, T imag)
{
  r = real;
  i = imag;
}

template <typename T>
VESSL_INLINE void transform::complex<T>::set_polar(T magnitude, T angle)
{
  r = magnitude * math::cos<T>(angle);
  i = magnitude * math::sin<T>(angle);
}

template<typename T>
VESSL_INLINE T transform::complex<T>::normalize()
{
  static constexpr T one = cast<T>(1);
  T m = math::sqrt(r*r + i*i);
  T d = m > math::epsilon<T>() ? one / m : one;
  r *= d;
  i *= d;
  return m;
}

template<typename T>
VESSL_INLINE T transform::complex<T>::magnitude() const
{
  return math::sqrt(r*r + i*i);
}

template <typename T>
VESSL_INLINE T transform::complex<T>::magnitude_sqr() const 
{ 
  return r*r + i*i; 
}

template <typename T>
phase_t transform::complex<T>::phase() const
{
  T phase_rad = math::atan2<T>(i, r);
  return cast<phase_t>(phase_rad / math::two_pi<T>());
}

template <typename T>
VESSL_INLINE transform::fft<T>::fft() { static_assert("fft is not available for this type"); }

template <typename T>
VESSL_INLINE transform::fft<T>::fft(size_t size) { static_assert("fft is not available for this type"); }

template <typename T>
VESSL_INLINE void transform::fft<T>::initialize(size_t size) {}

template <typename T>
VESSL_INLINE void transform::fft<T>::size() const {}

template <typename T>
VESSL_INLINE void transform::fft<T>::forward(array<sample_t> input, array<complex_t> output) {}

template <typename T> VESSL_INLINE void transform::fft<T>::inverse(array<complex_t> input, array<sample_t> output) {}
} // namespace vessl

#ifdef ARM_CORTEX
#include <arm_math.h>
#include "vessl.h"

namespace vessl
{
namespace transform
{
template<>
class fft<float32_t>
{
  arm_rfft_fast_instance_f32 rfft_instance_ = {};
  
public:
  using sample_t = float32_t;
  using complex_t = complex<sample_t>;
    
  VESSL_INLINE fft() = default;
  VESSL_INLINE explicit fft(size_t size)
  {
    initialize(size);
  }
    
  VESSL_INLINE void initialize(size_t size)
  {
    VASSERT(size==32 || size ==64 || size==128 || size==256 || size==512 || size==1024 || size==2048 || size==4096, "Unsupported FFT size");
    arm_rfft_fast_init_f32(&rfft_instance_, size);
  }
  
  VESSL_INLINE size_t size() const
  {
    return rfft_instance_.fftLenRFFT;
  }
  
  void forward(array<sample_t> input, array<complex_t> output)
  {
    VASSERT(input.size() >= size(), "Input array too small");
    VASSERT(output.size() >= size()/2, "Output array too small");
    arm_rfft_fast_f32(&rfft_instance_, input.data(), reinterpret_cast<sample_t *>(output.data()), 0);
  }
  
  void inverse(array<complex_t> input, array<sample_t> output)
  {
    VASSERT(input.size() >= size()/2, "Input array too small");
    VASSERT(output.size() >= size(), "Output array too small");
    arm_rfft_fast_f32(&rfft_instance_, reinterpret_cast<sample_t *>(input.data()), output.data(), 1);
  }
};
}
}

#endif