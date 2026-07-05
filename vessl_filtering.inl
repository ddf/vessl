#pragma once

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
