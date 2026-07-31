#pragma once

namespace vessl
{
template <typename T>
void transform::complex<T>::scale(T scalar)
{
  samples[0] *= scalar;
  samples[1] *= scalar;
}

template <typename T>
VESSL_INLINE void transform::complex<T>::set_complex(T real, T imag)
{
  samples[0] = real;
  samples[1] = imag;
}

template <typename T>
VESSL_INLINE void transform::complex<T>::set_polar(T magnitude, T angle)
{
  samples[0] = magnitude * math::cos<T>(angle);
  samples[1] = magnitude * math::sin<T>(angle);
}

template <typename T>
VESSL_INLINE T transform::complex<T>::magnitude() const
{
  return math::sqrt(samples[0]*samples[0] + samples[1]*samples[1]);
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

template <typename T>
VESSL_INLINE void transform::fft<T>::inverse(array<complex_t> input, array<sample_t> output) {}
} // namespace vessl

#ifdef ARM_CORTEX
#include <arm_math.h>

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