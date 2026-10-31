#pragma once

#include "vessl.h"

template <typename T>
VESSL_INLINE void vessl::filtering::dc_block<T>::process(const T *source, T *dest, size_t block_size, const args &args)
{
  // gives us about 0.995 for 44100, which is a pretty good R according to the article above.
  analog_t r = (args.sr - 200.0f) / args.sr;
  while (block_size--)
  {
    T x = *source++;
    y1 = x - x1 + r*y1;
    x1 = x;
    *dest++ = y1;
  }
}

template <typename T, vessl::filtering::moog::type F>
VESSL_INLINE void vessl::filtering::moog::flt1<T, F>::process(const T* source, T* dest, size_t block_size, const vessl::filtering::args& args)
{
  analog_t fn = args.hz / (args.sr*0.5f);
  analog_t rz = vessl::math::constrain(args.q, 0.f, 1.f);

  analog_t q = 1.f - fn;
  analog_t p = fn + 0.8f * fn * q;
  analog_t f = p + p - 1.f;
  q = rz * (1.f + 0.5f * q * (1.f - q + 5.6f * q * q));

  analog_t t1, t2;
  while(block_size--)
  {
    analog_t in = vessl::cast<analog_t>(*source++);
    in -= q * b4; // feedback
    
    t1 = b1;
    b1 = (in + b0)*p - b1*f;

    t2 = b2;
    b2 = (b1 + t1)*p - b2*f;

    t1 = b3;
    b3 = (b2 + t2)*p - b3*f;

    b4 = (b3 + t1)*p - b4*f;
    b4 = b4 - b4*b4*b4*b4*0.166667f; // clipping
    
    // inelegantly squash denormals
    if (vessl::math::is_nan(b4))
    {
      b1 = b2 = b3 = b4 = 0;
    }

    b0 = in;

    if constexpr(F == vessl::filtering::moog::type::lp)
    {
      *dest++ = vessl::cast<T>(b4);
    }

    if constexpr(F == vessl::filtering::moog::type::bp)
    {
      *dest++ = vessl::cast<T>(3.f * (b3-b4)); 
    }

    if constexpr(F == vessl::filtering::moog::type::hp)
    {
      *dest++ = vessl::cast<T>(in - b4);
    }
  }
}

template <typename T>
VESSL_INLINE void vessl::filtering::moog::flt2<T>::process(const T* source, T* dest, size_t block_size, const vessl::filtering::args& args)
{
  analog_t fc = args.hz / (args.sr*0.5f);
  analog_t rz = vessl::math::constrain(args.q, 0.f, 4.f);

  analog_t f = fc * 1.16f;
  analog_t fb = rz * (1.f - 0.15f * f * f);
  while(block_size--)
  {
    analog_t input = vessl::cast<analog_t>(*source++);
    input -= out4 * fb;
    input *= 0.35013 * (f*f)*(f*f);
    // the 0.3 scaling of the input is what keeps resonance under control when in the [0,4] range.
    out1 = input + 0.3 * in1 + (1 - f) * out1; // Pole 1
    in1  = input;
    out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
    in2  = out1;
    out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
    in3  = out2;
    out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
    in4  = out3;
    *dest++ = vessl::cast<T>(out4);
  }
}

template <typename T>
VESSL_INLINE vessl::filtering::data<T>::data(coeff_t *coeff_data, size_t coeff_size, T *state_data, size_t state_size)
: coeff(coeff_data)
, state(state_data)
{
  array<T> carr(coeff_data, coeff_size);
  carr.fill(0);
  
  array<T> sarr(state_data, state_size);
  sarr.fill(T(0));
}

namespace vessl
{
template<size_t Stages>
template<typename T, class CoGen>
VESSL_INLINE void filtering::biquad<Stages>::df2t<T, CoGen>::process(const T* source, T* dest, size_t block_size, const args& args)
{
  // update our coefficients
  using cascade_t = cascade<T,2>;
  using coeff_t = typename cascade_t::coeff_t;
  
  cg(cascade_t::coeff, args.omega(), args.q, args.g);
    
  // run the filter
  const T* input = source;
  typename array<coeff_t>::reader cor(cascade_t::coeff, coeff_count);
  for (size_t s = 0; s < Stages; s++)
  {
    analog_t b0 = cor.read(); analog_t b1 = cor.read(); analog_t b2 = cor.read();
    analog_t a1 = cor.read();  analog_t a2 = cor.read();
    T* st = cascade_t::state + 2*s;
    T d1 = st[0]; T d2 = st[1];
    typename array<T>::reader r(input, block_size);
    typename array<T>::writer w(dest, block_size);
    while (r)
    {
      T xn = r.read();
      T yn = b0 * xn + d1;
      d1 = b1 * xn + a1 * yn + d2;
      d2 = b2 * xn + a2 * yn;
      w << yn;
    }
    st[0] = d1;
    st[1] = d2;
    input = dest;
  }
}

template<size_t Stages>
template<typename T>
VESSL_INLINE void filtering::biquad<Stages>::copy(T* coeff) 
{
  if (Stages > 1)
  {
    array<T> src(coeff, coeff_num);
    for (size_t i = 1; i < Stages; i++)
    {
      array<T> dst(coeff + coeff_num*i, coeff_num);
      src.copy_to(dst);
    }
  }
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::lpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
{
  analog_t K = math::tan(omega);
  analog_t A = 1 / (1 + K / q + K * K);
  coeff[0] = K * K * A;
  coeff[1] = 2 * coeff[0];
  coeff[2] = coeff[0];
  coeff[3] = - 2 * (K * K - 1) * A;
  coeff[4] = - (1 - K / q + K * K) * A;
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::hpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
{
  analog_t K = math::tan(omega);
  analog_t A = 1 / (1 + K / q + K * K);
  coeff[0] = 1 * A;
  coeff[1] = -2 * coeff[0];
  coeff[2] = coeff[0];
  coeff[3] = -2 * (K * K - 1) * A;
  coeff[4] = -(1 - K / q + K * K) * A;
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::bpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
{
  analog_t K = math::tan(omega);
  analog_t A = 1 / (1 + K / q + K * K);
  coeff[0] = K / q * A;
  coeff[1] = 0;
  coeff[2] = -coeff[0];
  coeff[3] = -2 * (K * K - 1) * A;
  coeff[4] = -(1 - K / q + K * K) * A;
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::ntcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
{
  analog_t K = math::tan(omega);
  analog_t A = 1 / (1 + K / q + K * K);
  coeff[0] = (1 + K * K) * A;
  coeff[1] = 2 * (K * K - 1) * A;
  coeff[2] = coeff[0];
  coeff[3] = -coeff[1];
  coeff[4] = - (1 - K / q + K * K) * A;
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::pkcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t g) const
{
  analog_t K = math::tan(omega);
  analog_t V = math::exp10(math::abs(g.to_decibels())/20);
  analog_t A;
  if (g)
  {
    A = 1 / (1 + K / q + K * K);
    coeff[0] = (1 + V/q * K + K * K) * A;
    coeff[1] = 2 * (K * K - 1) * A;
    coeff[2] = (1 - V/q * K + K * K) * A;
    coeff[3] = -coeff[1];
    coeff[4] = - (1 - K / q + K * K) * A;
  }
  else
  {
    A = 1 / (1 + V/q * K + K * K);
    coeff[0] = (1 + K/q + K * K) * A;
    coeff[1] = 2 * (K * K - 1) * A;
    coeff[2] = (1 - K/q + K * K) * A;
    coeff[3] = -coeff[1];
    coeff[4] = - (1 - V/q * K + K * K) * A;
  }
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::lscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const
{
  analog_t K = math::tan(omega);
  analog_t V = math::exp10(math::abs(g.to_decibels())/20);
  analog_t A;
  if (g)
  {
    A = 1 / (1 + math::sqrt2<analog_t>() * K + K * K);
    coeff[0] = (1 + math::sqrt(2*V) * K + V * K * K) * A;
    coeff[1] = 2 * (V * K * K - 1) * A;
    coeff[2] = (1 - math::sqrt(2*V) * K + V * K * K) * A;
    coeff[3] = -2 * (K * K - 1) * A;
    coeff[4] = -(1 - math::sqrt2<analog_t>()*K + K * K) * A;
  }
  else
  {
    A = 1 / (1 + math::sqrt(2*V) * K + V * K * K);
    coeff[0] = (1 + math::sqrt2<analog_t>()*K + K * K) * A;
    coeff[1] = 2 * (K * K - 1) * A;
    coeff[2] = (1 - math::sqrt2<analog_t>()*K + K * K) * A;
    coeff[3] = -2 * (V * K * K - 1) * A;
    coeff[4] = -(1 - math::sqrt(2*V) * K + V * K * K) * A;
  }
  copy<analog_t>(coeff);
}

template<size_t Stages>
VESSL_INLINE void filtering::biquad<Stages>::hscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const
{
  analog_t K = math::tan(omega);
  analog_t V = math::exp10(math::abs(g.to_decibels())/20);
  analog_t A;
  if (g)
  {
    A = 1 / (1 + math::sqrt2<analog_t>() * K + K * K);
    coeff[0] = (V + math::sqrt(2*V) * K + K * K) * A;
    coeff[1] = 2 * (K * K - V) * A;
    coeff[2] = (V - math::sqrt(2*V) * K + K * K) * A;
    coeff[3] = -2 * (K * K - 1) * A;
    coeff[4] = -(1 - math::sqrt2<analog_t>()*K + K * K) * A;
  }
  else
  {
    A = 1 / (V + math::sqrt(2*V) * K + K * K);
    coeff[0] = (1 + math::sqrt2<analog_t>()*K + K * K) * A;
    coeff[1] = 2 * (K * K - 1) * A;
    coeff[2] = (1 - math::sqrt2<analog_t>()*K + K * K) * A;
    coeff[3] = -2 * (K * K - V) * A;
    coeff[4] = -(V - math::sqrt(2*V) * K + K * K) * A;
  }
  copy<analog_t>(coeff);
}
} // namespace vessl
