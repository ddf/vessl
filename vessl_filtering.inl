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

namespace vessl::filtering::moog
{
// table to fixup Q in order to remain constant for various pole frequencies, from Tim Stilson's code @ CCRMA (also in CLM distribution)
static const float gain_table[199] = { 
  0.999969, 0.990082, 0.980347, 0.970764, 0.961304, 0.951996, 0.94281, 0.933777, 0.924866, 0.916077, 
  0.90741, 0.898865, 0.890442, 0.882141 , 0.873962, 0.865906, 0.857941, 0.850067, 0.842346, 0.834686, 
  0.827148, 0.819733, 0.812378, 0.805145, 0.798004, 0.790955, 0.783997, 0.77713, 0.770355, 0.763672, 
  0.75708 , 0.75058, 0.744141, 0.737793, 0.731537, 0.725342, 0.719238, 0.713196, 0.707245, 0.701355, 
  0.695557, 0.689819, 0.684174, 0.678558, 0.673035, 0.667572, 0.66217, 0.65686, 0.651581, 0.646393, 
  0.641235, 0.636169, 0.631134, 0.62619, 0.621277, 0.616425, 0.611633, 0.606903, 0.602234, 0.597626, 
  0.593048, 0.588531, 0.584045, 0.579651, 0.575287 , 0.570953, 0.566681, 0.562469, 0.558289, 0.554169, 
  0.550079, 0.546051, 0.542053, 0.538116, 0.53421, 0.530334, 0.52652, 0.522736, 0.518982, 0.515289, 
  0.511627, 0.507996 , 0.504425, 0.500885, 0.497375, 0.493896, 0.490448, 0.487061, 0.483704, 0.480377, 
  0.477081, 0.473816, 0.470581, 0.467377, 0.464203, 0.46109, 0.457977, 0.454926, 0.451874, 0.448883, 
  0.445892, 0.442932, 0.440033, 0.437134, 0.434265, 0.431427, 0.428619, 0.425842, 0.423096, 0.42038, 
  0.417664, 0.415009, 0.412354, 0.409729, 0.407135, 0.404572, 0.402008, 0.399506, 0.397003, 0.394501, 
  0.392059, 0.389618, 0.387207, 0.384827, 0.382477, 0.380127, 0.377808, 0.375488, 0.37323, 0.370972, 
  0.368713, 0.366516, 0.364319, 0.362122, 0.359985, 0.357849, 0.355713, 0.353607, 0.351532, 0.349457, 
  0.347412, 0.345398, 0.343384, 0.34137, 0.339417, 0.337463, 0.33551, 0.333588, 0.331665, 0.329773, 
  0.327911, 0.32605, 0.324188, 0.322357, 0.320557, 0.318756, 0.316986, 0.315216, 0.313446, 0.311707, 
  0.309998, 0.308289, 0.30658, 0.304901, 0.303223, 0.301575, 0.299927, 0.298309, 0.296692, 0.295074, 
  0.293488, 0.291931, 0.290375, 0.288818, 0.287262, 0.285736, 0.284241, 0.282715, 0.28125, 0.279755, 
  0.27829, 0.276825, 0.275391, 0.273956, 0.272552, 0.271118, 0.269745, 0.268341, 0.266968, 0.265594, 
  0.264252, 0.262909, 0.261566, 0.260223, 0.258911, 0.257599, 0.256317, 0.255035, 0.25375 
};

template <typename T, type F>
VESSL_INLINE void flt<T, F>::process(const T* source, T* dest, size_t block_size, const args& args)
{
  /** @todo improve frequency tracking
   *  tracking increasingly falls short below ~1.5k and starts to overshoot above ~2.1k
   */
  analog_t fc = args.fc();
  analog_t rz = vessl::math::constrain(args.q, 0.f, 1.f);
  analog_t x2 = fc*fc;
  analog_t x3 = fc*x2;
  analog_t p  = -0.69346 * x3 - 0.59515 * x2 + 3.2937 * fc - 1.0072; //cubic fit by DFL, not 100% accurate but better than nothing...
  // when sending analog noise thru this filter, we found that it needed to be attenuated by 12db to prevent distortion
  // so we use that as a baseline and gain from the args can adjust if desired.
  analog_t g  = 0.25f * args.g.to_scale();

  analog_t ix, ixfrac;
  int ixint;
  ix = p * 99;
  ixint = vessl::math::floor(ix);
  ixfrac = ix - ixint;
  analog_t Q;
  vessl::sample::crossfade( gain_table[ ixint + 99 ], gain_table[ ixint + 100 ], ixfrac, &Q );
  Q *= args.q;

  // note: tried using softlimit instead of saturate, which reduces distortion when filtering noise,
  // but nerfs the ping response, which we'd like to keep.
  while(block_size--)
  {
    analog_t input = vessl::cast<analog_t>(*source++)*g;
    output = 0.25 * ( input - output ); //negative feedback

    analog_t temp = b0;
    output = vessl::sample::saturate( output + p * (output - temp));
    b0 = output;
    output = vessl::sample::saturate( output + temp );

    temp = b1;
    output = vessl::sample::saturate( output + p * (output - temp));
    b1 = output;
    output = vessl::sample::saturate( output + temp );

    temp = b2;
    output = vessl::sample::saturate( output + p * (output - temp));
    b2 = output;
    output = vessl::sample::saturate( output + temp );

    temp = b3;
    output = vessl::sample::saturate( output + p * (output - temp));
    b3 = output;
    output = vessl::sample::saturate( output + temp );

    if constexpr(F == type::lp)
    {
      *dest++ = vessl::cast<T>(output);
    }

    if constexpr(F == type::bp)
    {
      *dest++ = vessl::cast<T>(3 * (b3 - output)); //got this one from paul kellet
    }

    if constexpr(F == type::hp)
    {
      *dest++ = vessl::cast<T>(input - output);
    }

    output *= Q;  //scale the feedback
  }

}

template <typename T, type F>
VESSL_INLINE void v1::flt<T, F>::process(const T* source, T* dest, size_t block_size, const args& args)
{
  /** fc scaling is our addition, to improve the frequency tracking of the filter.
   */
  analog_t c = args.fc() * 0.899f;
  analog_t r = vessl::math::constrain(args.q, 0.f, 1.f);
  analog_t g = args.g.to_scale();

  analog_t d = 1.f - c;
  // original
  analog_t p = (c + 0.8f*c*d);
  // our modified version - not sure it's better?
  // see: https://www.desmos.com/calculator/nccplgcgrh
  //analog_t p = c*(d + 0.8f*d*d + 0.25f);
  analog_t f = p + p - 1.f;
  analog_t q = r * (1.f + 0.5f * d * (c + 5.6f*(d*d)));

  analog_t t1, t2;
  while(block_size--)
  {
    analog_t in = vessl::cast<analog_t>(*source++)*g;
    in -= q * b4; // feedback
    
    t1 = b1; b1 = (in + b0)*p - b1*f;
    t2 = b2; b2 = (b1 + t1)*p - b2*f;
    t1 = b3; b3 = (b2 + t2)*p - b3*f;
             b4 = (b3 + t1)*p - b4*f;
    b4 = b4 - b4*b4*b4*0.166667f; // clipping
    
    // inelegantly squash denormals
    if (vessl::math::is_nan(b4))
    {
      b1 = b2 = b3 = b4 = 0;
    }

    b0 = in;

    if constexpr(F == type::lp)
    {
      *dest++ = vessl::cast<T>(b4);
    }

    if constexpr(F == type::bp)
    {
      *dest++ = vessl::cast<T>(3.f * (b3-b4)); 
    }

    if constexpr(F == type::hp)
    {
      *dest++ = vessl::cast<T>(in - b4);
    }
  }
}

template <typename T, type F>
VESSL_INLINE void v2::flt<T, F>::process(const T* source, T* dest, size_t block_size, const args& args)
{
  // the scaling up of fc is our addition.
  // without it, this filter tracks far below the cutoff frequency,
  // especially compared with v1.
  analog_t fc = args.fc() * 2.62f;
  analog_t rz = vessl::math::constrain(args.q, 0.f, 1.f)*3.99f;
  analog_t gn = args.g.to_scale();

  analog_t f = fc * 1.16f;
  analog_t fb = rz * (1.f - 0.15f * f * f);
  while(block_size--)
  {
    analog_t input = vessl::cast<analog_t>(*source++)*gn;
    input -= out4 * fb; // need this for hp
    
    analog_t in0 = input * 0.35013 * (f*f)*(f*f);
    // the 0.3 scaling of the input is what keeps resonance under control when in the [0,4] range.
    out1 = in0 + 0.3 * in1 + (1.f - f) * out1; // Pole 1
    in1  = in0;
    out2 = out1 + 0.3 * in2 + (1.f - f) * out2;  // Pole 2
    in2  = out1;
    out3 = out2 + 0.3 * in3 + (1.f - f) * out3;  // Pole 3
    in3  = out2;
    out4 = out3 + 0.3 * in4 + (1.f - f) * out4;  // Pole 4
    in4  = out3;

    if constexpr(F == type::lp)
    {
      *dest++ = vessl::cast<T>(out4);
    }

    if constexpr(F == type::bp)
    {
      *dest++ = vessl::cast<T>(3.f * (out3 - out4)); 
    }

    if constexpr(F == type::hp)
    {
      *dest++ = vessl::cast<T>(input - out4);
    }
  }
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
