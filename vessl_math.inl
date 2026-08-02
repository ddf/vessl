#pragma once

#ifdef ARM_CORTEX
#include "vessl_arm_math.inl"
#endif

namespace vessl
{
namespace math
{

VESSL_INLINE analog_t decibels_to_scale(analog_t db)
{
  return exp10(db*cast<analog_t>(0.05));
}

VESSL_INLINE analog_t scale_to_decibels(analog_t scale)
{
  return log10(scale)*cast<analog_t>(20.0);
}

template<>
VESSL_INLINE analog_t sin<analog_t, phase_t>(phase_t z) 
{ 
  return sin<analog_t>(two_pi<analog_t>() * cast<analog_t>(z)); 
}

template<>
VESSL_INLINE analog_t cos<analog_t, phase_t>(phase_t z) 
{ 
  return cos<analog_t>(two_pi<analog_t>() * cast<analog_t>(z)); 
}

template <>
VESSL_INLINE analog_t mod<analog_t>(analog_t v, analog_t *i)
{
  return ::modff(v, i);
}

template <typename T>
VESSL_INLINE T abs(const T &val)
{
  return val < 0 ? -val : +val;
}

template<>
VESSL_INLINE float abs<float>(const float &val)
{
  return fabsf(val);
}

template<>
VESSL_INLINE double abs<double>(const double &val)
{
  return fabs(val);
}

template <typename T>
constexpr VESSL_INLINE T epsilon()
{
  return T(0);
}

template<>
constexpr VESSL_INLINE float epsilon<float>()
{
  return FLT_EPSILON;
}

template<>
constexpr VESSL_INLINE double epsilon<double>()
{
  return DBL_EPSILON;
}
  
template<typename T, size_t N>
VESSL_INLINE sample::frame<T, N> round(sample::frame<T, N> x)
{
  sample::frame<T, N> result;
  for (size_t i = 0; i < N; i++)
  {
    result.samples[i] = round(x.samples[i]);
  }
  return result;
}

template<>
VESSL_INLINE analog_t xore(const analog_t& a, const analog_t& b)
{
  return xore(cast<digital_t>(a), cast<digital_t>(b));
}

namespace random
{
VESSL_INLINE uint32_t u32()
{
  ru32_seed ^= ru32_seed << 13; ru32_seed ^= ru32_seed >> 17; ru32_seed ^= ru32_seed << 5;
  return ru32_seed;
}

template<typename T>
VESSL_INLINE T range(T low, T high)
{
  static constexpr analog_t scale = 1/4294967296.0;
  analog_t r = static_cast<analog_t>(u32()) * scale;
  return low + r*(high-low);
}
} // namespace random

namespace easing
{

template<typename D>
VESSL_INLINE D linear::operator()(D t) const
{
  return t;
}

VESSL_INLINE analog_t smoothstep::operator()(analog_t t) const
{
  return t * t * (3.0f - 2.0f * t);
}

VESSL_INLINE analog_t quart::in::operator()(analog_t t) const
{
  return t*t*t*t;
}

VESSL_INLINE analog_t quart::out::operator()(analog_t t) const
{
  return 1.f - pow<analog_t>(1.f - t, 4.f);
}

VESSL_INLINE analog_t quart::in_out::operator()(analog_t t) const
{
  return t < 0.5f ? 8.f * t * t * t * t : 1.f - pow<analog_t>(-2.f * t + 2.f, 4.f) * 0.5f;
}

VESSL_INLINE analog_t expo::in::operator()(analog_t t) const
{
  return t <= epsilon<analog_t>() ? 0 : pow<analog_t>(2, 10*t-10);
}

VESSL_INLINE analog_t expo::out::operator()(analog_t t) const
{
  return t >= 1.0 - epsilon<analog_t>() ? 1.0 : 1.0 - pow<analog_t>(2, -10*t);
}

VESSL_INLINE analog_t expo::in_out::operator()(analog_t t) const
{
  return t <= math::epsilon<analog_t>() ? 0.0
           : t >= 1.0 - math::epsilon<analog_t>() ? 1.0  // NOLINT(cppcoreguidelines-narrowing-conversions)
           : t < 0.5 ? math::pow<analog_t>(2.0, 20*t-10)*0.5  // NOLINT(clang-diagnostic-implicit-float-conversion)
           : 2.0 - math::pow<analog_t>(2.0, -20*t+10)*0.5;  // NOLINT(clang-diagnostic-implicit-float-conversion)
}

VESSL_INLINE analog_t quad::in::operator()(analog_t t) const
{
  return t * t;
}

VESSL_INLINE analog_t quad::out::operator()(analog_t t) const
{
  return 1.0 - (1.0 - t) * (1.0 - t);  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions)
}

VESSL_INLINE analog_t quad::in_out::operator()(analog_t t) const
{
  return t < 0.5 ? 2*t*t : 1.0 - pow<analog_t>(-2*t+2, 2) * 0.5;
}

VESSL_INLINE analog_t quad::out_in::operator()(analog_t t) const
{
  static out qo; 
  return t < 0.5 ? qo(2*t) * 0.5f : 1.0f - qo(2*t) * 0.5f;
}

template <typename T>
VESSL_INLINE T smooth(T value, T target, analog_t degree)
{
  return value*degree + (1.0 - degree)*target;
}

template <typename T>
smoother<T>::smoother(analog_t smoothing_degree, T initial_value) 
: value(initial_value), degree(smoothing_degree)
{
}

template <typename T>
VESSL_INLINE smoother<T>& smoother<T>::operator=(const T &v)
{
  analog_t d = constrain<analog_t>(degree, 0.0, 1.0);
  value = smooth(value, v, d);
  return *this;
}

template<>
VESSL_INLINE digital_t smooth<digital_t>(digital_t value, digital_t target, analog_t degree)
{
  return (value*degree + target)/(degree+1); // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions)
}
} // namespace easing

template<typename D>
VESSL_INLINE D clamp_delta(D t)
{
  return t < 0 ? 0 : t > 1 ? 1 : t;
}

template<>
VESSL_INLINE phase_t clamp_delta<phase_t>(phase_t t)
{
  return t;
}

template<>
VESSL_INLINE q31 clamp_delta<q31>(q31 t)
{
  return t < q31::mid() ? q31::mid() : t > q31::max() ? q31::max() : t;
}

template<typename E, typename T, typename D>
VESSL_INLINE T interp(T begin, T end, D t)
{
  static E ease;
  t = clamp_delta(t);
  return (end-begin) * ease(t) + begin;
}

template<typename T, typename D>
struct lerper
{
  VESSL_INLINE T operator()(T begin, T end, D t)
  {
    return interp<easing::linear>(begin, end, t);
  }
};

template<typename T>
struct lerper<T,phase_t>
{
  VESSL_INLINE T operator()(T begin, T end, phase_t t)
  {
    return t == phase_zero ? begin 
         : t == phase_360 ? end 
         : begin < end ? begin + (end-begin)*t/phase_360
         : begin - (begin-end)*t/phase_360; 
  }
};

template <typename T>
T atan2(T y, T x)
{
  return cast<T>(::atan2(cast<double>(y), cast<double>(x)));
}

template <typename T, typename D>
VESSL_INLINE T lerp(T begin, T end, D t)
{
  static lerper<T,D> func;
  return func(begin, end, t);
}

template<typename T>
VESSL_INLINE T wrap(T val, T low, T high)
{
  // @todo probably a way to do this without while loops.
  T diff = high - low;
  while (val < low) { val += diff; }
  while (val > high) { val -= diff; }
  return val; 
}

template<>
VESSL_INLINE phase_t wrap01<phase_t>(phase_t val)
{
  return val;
}
      
template<>
VESSL_INLINE analog_t wrap01<analog_t>(analog_t v)
{
  analog_t i; 
  analog_t f = mod(v, &i);
  return f < 0 ? f + 1.0f : f;
}

} // namespace math
} // namespace vessl

namespace vessl
{

template <typename T>
VESSL_INLINE transform33<T>::transform33(): matrix<T>(data_, 3, 3)
{
  set_identity();
}

template <typename T>
transform33<T>::transform33(const transform33 &other): matrix<T>(data_, 3, 3)
{
  memcpy(data_, other.data_, sizeof(T)*3*3);
}

template <typename T>
VESSL_INLINE void transform33<T>::set_identity()
{
  matrix<T>::clear();
  for (size_t i = 0; i < 3; i++) {
    matrix<T>::set(i,i, T(1LL));
  }
}

template <typename T>
VESSL_INLINE void transform33<T>::set_euler(phase_t pitch, phase_t yaw, phase_t roll)
{
  T cosa = math::cos<T>(roll);
  T sina = math::sin<T>(roll);

  T cosb = math::cos<T>(yaw);
  T sinb = math::sin<T>(yaw);

  T cosc = math::cos<T>(pitch);
  T sinc = math::sin<T>(pitch);
  
  using m = matrix<T>;

  // row*3 + col
  m::set(0, 0, cosa * cosb);
  m::set(0, 1, cosa * sinb*sinc - sina*cosc);
  m::set(0, 2, cosa * sinb*cosc + sina*sinc);

  m::set(1, 0, sina * cosb);
  m::set(1, 1, sina * sinb*sinc + cosa*cosc);
  m::set(1, 2, sina * sinb*cosc - cosa*sinc);

  m::set(2, 0, -sinb);
  m::set(2, 1, cosb * sinc);
  m::set(2, 2, cosb * cosc);
}

template <typename T>
VESSL_INLINE void transform33<T>::set_euler_radians(analog_t pitch_radians, analog_t yaw_radians, analog_t roll_radians)
{
  return set_euler(cast<phase_t>(pitch_radians / math::two_pi<analog_t>()), 
                   cast<phase_t>(yaw_radians / math::two_pi<analog_t>()), 
                   cast<phase_t>(roll_radians / math::two_pi<analog_t>()));
}

template <typename T>
VESSL_INLINE sample::frame<T, 3> transform33<T>::multiply(const sample::frame<T, 3> &input)
{
  using m = matrix<T>;
  sample::frame<T,3> output;
  T* out = output.samples;
  const T* in = input.samples;

  out[0] = m::get(0,0) * in[0] + m::get(0,1) * in[1] + m::get(0,2) * in[2];
  out[1] = m::get(1,0) * in[0] + m::get(1,1) * in[1] + m::get(1,2) * in[2];
  out[2] = m::get(2,0) * in[0] + m::get(2,1) * in[1] + m::get(2,2) * in[2];

  // this might be faster?
  //mtrx.multiply(input.toMatrix(), output.toMatrix());

  return output;
}
}