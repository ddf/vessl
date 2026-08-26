#pragma once

#include "vessl_units.h"

namespace vessl
{
namespace generators
{
template<typename T, typename N>
VESSL_INLINE T noise<T, N>::generate()
{
  step_ += dt_ * params_.rate.value;
  if (step_ >= 1)
  {
    value_ = next_;
    next_ = noise_source();
    step_ = math::wrap01(step_);
  }
  return value_;
}

template<typename T, typename N>
template<typename E>
VESSL_INLINE T noise<T, N>::generate()
{
  T s = generate();
  return math::interp<E>(s, next_, step_);
}

template<typename T>
VESSL_INLINE void ramp<T>::trigger()
{
  if (duration() > math::epsilon<analog_t>())
  {
    t_ = 0;
    params_.eor.value = false;
  }
  else
  {
    t_ = 1;
    params_.eor.value = true;
  }
}

template<typename T>
VESSL_INLINE T ramp<T>::generate()
{
  analog_t lt = t_;
  if (is_active())
  {
    analog_t dinv = 1.f / duration();
    t_ += dt_ * dinv;
    if (t_ >= 1)
    {
      params_.eor.value = true;
    }
  }
  return math::lerp(params_.from.value, params_.to.value, lt);
}

template<typename T>
template<typename E>
VESSL_INLINE T envelope<T>::stage::step()
{
  analog_t s = dt_ / math::max<analog_t>(cast<analog_t>(duration()), dt_);
  time_ += s;
  analog_t t = math::constrain<analog_t>(time_, 0.0, 1.0);
  params_.output.value = math::interp<E>(begin_, params_.target.value, t);
  if (time_ >= 1)
  {
    params_.active.value = false;
    params_.eos.value = true;
  }
  return params_.output.value;
}

template<typename T>
VESSL_INLINE void envelope<T>::trigger()
{
  for (stage& s : stages_)
  {
    s.reset();
  }
  final_.reset();
  params_.eoc.value = false;
  stage_idx_ = 0;
  stages_[0].start(0);
}

template<typename T>
template<typename E>
VESSL_INLINE T envelope<T>::generate()
{
  T value = current_stage().template generate<E>();
  if (stage_idx_ == stages_.size() && final_.eos())
  {
    params_.eoc.value = true;
  }
  else if (should_advance(stage_idx_))
  {
    get_stage(++stage_idx_).start(value);
  }
  return value;
}

template<typename T>
VESSL_INLINE void envelope<T>::set_sample_rate(analog_t sample_rate)
{
  for (stage& s : stages_)
  {
    s.set_sample_rate(sample_rate);
  }
  final_.set_sample_rate(sample_rate);
}

template<typename T>
VESSL_INLINE void asr<T>::gate(T value)
{
  binary_t value_on = value > trig_threshold_;
  T attack_target = attack().target().template read<T>();
  if (!gate_on_ && value_on)
  {
    ad<T>::trigger();
    gate_on_ = true;
    // reset the attackTarget because where the gate value ends up might be less than the previous attackTarget
    attack_target = 0;
  }
  else if (gate_on_ && !value_on)
  {
    gate_on_ = false;
    // if we're mid-attack, start the release
    if (attack().active())
    {
      envelope<T>::start_stage(1, attack().value().template read<T>());
    }
  }
  attack().target() = math::max(value, attack_target);
}

template<class W>
VESSL_INLINE typename W::sample_t oscil<W>::generate()
{
  typename W::sample_t val = waveform.evaluate(phase_ + params_.pm.value);
  analog_t f = params_.fhz.value * math::exp2(params_.fm_exp.value) + params_.fm_lin.value;
  phase_ += static_cast<phase_t>(dt_ * f);
  return val;
}

template <class W>
VESSL_INLINE void oscil<W>::generate(sink<typename W::sample_t>& dest)
{
  analog_t freq = params_.fHz.value * math::exp2(params_.fm_exp.value) + params_.fm_lin.value;
  phase_t step = static_cast<phase_t>(dt_ * freq);
  while(!dest.isFull())
  {
    dest << waveform.evaluate(phase_ + params_.pm.value);
    phase_ += step;
  }
}

template <typename T>
clock<T>::clock(analog_t sample_rate, period_t sample_period_min, period_t sample_period_max, analog_t bpm)
  : clockable(sample_rate, sample_period_min, sample_period_max, bpm)
  , dt_(cast<phase_t>(1.0f/sample_rate)), phase_(phase_zero)
  , pulse_(phase_90)
{
    
}

template <typename T>
VESSL_INLINE parameter clock<T>::tempo() const
{
  return parameter("tempo",'t', &tempo_);
}

template <typename T>
VESSL_INLINE T clock<T>::generate()
{
  tick();
  T val = pulse_.evaluate(phase_);
  analog_t freq = tempo_.to_frequency(sample_rate_);
  phase_ += static_cast<phase_t>(dt_ * freq);
  return val;
}

template <typename T>
VESSL_INLINE parameter clock<T>::element_at(size_t index) const
{
  return tempo();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// spectral
//
template <typename T, size_t SpectrumSize, size_t Overlap>
constexpr spectral<T, SpectrumSize, Overlap>::frequency_band::frequency_band(T magnitude, phase_t phase)
: magnitude_(magnitude)
{
  complex_.set_polar(1.f, phase);
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE typename spectral<T, SpectrumSize, Overlap>::frequency_band& 
spectral<T, SpectrumSize, Overlap>::frequency_band::operator=(const frequency_band &other)
{
  if (this == &other)
  {
    return *this;
  }
  
  complex_ = other.complex_;
  magnitude_ = other.magnitude_;
  return *this;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE typename spectral<T, SpectrumSize, Overlap>::complex_t 
spectral<T, SpectrumSize, Overlap>::frequency_band::to_complex() const
{
  complex_t ret(complex_);
  ret.scale(magnitude_);
  return ret;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::set_complex(const complex_t& from_complex)
{
  complex_ = from_complex;
  magnitude_ = complex_.normalize();
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::set_polar(T magnitude, phase_t phase)
{
  magnitude_ = magnitude;
  complex_.set_polar(1.f, phase);
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::set_magnitude(T magnitude)
{
  magnitude_ = magnitude;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::add(const frequency_band &other)
{
  complex_t lhs = to_complex();
  complex_t rhs = other.to_complex();
  complex_.r = lhs.r + rhs.r;
  complex_.i = lhs.r + rhs.r;
  magnitude_ = complex_.normalize();
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::subtract(const frequency_band &other)
{
  complex_t lhs = to_complex();
  complex_t rhs = other.to_complex();
  complex_.r = lhs.r - rhs.r;
  complex_.i = lhs.r - rhs.r;
  magnitude_ = complex_.normalize();
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::frequency_band::blend(const frequency_band &other, analog_t amt)
{
  complex_.r = vessl::math::lerp(complex_.r, other.complex_.r, amt);
  complex_.i = vessl::math::lerp(complex_.i, other.complex_.i, amt);
  magnitude_ = vessl::math::lerp(magnitude_, other.magnitude_, amt);
}

template <typename T, size_t SpectrumSize, size_t Overlap>
spectral<T, SpectrumSize, Overlap>::spectral(data &data, analog_t sample_rate)
  : unit_generator<T>()
  , fft_(SpectrumSize)
  , bands_(data.bands.data(), data.bands.size())
  , spectrum_(data.spectrum.data(), data.spectrum.size())
  , window_(data.window.data(), data.window.size())
  , signal_(data.signal.data(), data.signal.size())
  , output_(data.output.data(), data.output.size())
  , read_idx_(0)
  , overlap_count_(block_size) // will trigger generation immediately the first time generate() is called
  , bin_spacing_(sample_rate/SpectrumSize)
  , phase_flip_(false)
{
  constexpr size_t band_count = SpectrumSize/2;
  VASSERT(data.bands.size() == band_count, "Invalid frequency bands size");
  VASSERT(data.spectrum.size() == band_count, "Invalid spectrum size");
  VASSERT(data.window.size() == SpectrumSize, "Invalid window size");
  VASSERT(data.signal.size() == SpectrumSize, "Invalid signal size");
  VASSERT(data.output.size() == SpectrumSize, "Invalid output size");

  bands_[0].set_complex({0,0});
  for (size_t i = 1; i < band_count; ++i)
  {
    frequency_band& band = bands_[i];
    band.set_polar(0, math::random::u32()/2);
  }

  output_.fill(0);
}

// for reference, see: https://dsp.stackexchange.com/questions/59519/definition-of-the-dft-fft-bin-size
template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE analog_t spectral<T, SpectrumSize, Overlap>::get_band_frequency(size_t index) const
{
  const analog_t k = static_cast<analog_t>(index);
  return k*bin_spacing_;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE size_t spectral<T, SpectrumSize, Overlap>::get_band_index(analog_t frequency) const
{
  const size_t k = static_cast<size_t>(math::round(frequency/bin_spacing_));
  return k;
}

/** 
spectrum:
  - SpectrumSize/2
  - can be shared between instances if processing is single-threaded
[--------------------------------]

signal generated by inverse
  - SpectrumSize
  - example overlap generation points marked with letters
  - can be shared if processing is single-threaded)
[a---b---c---d---e---f---g---h---a---b---c---d---e---f---g---h---]

mix buffer:
  - SpectrumSize
  - can be circular, but need to replace read samples with 0 as we go.
  - each overlap adds into the mix buffer starting at the next read position
[a---------------------------------------------------------------]
[----b-----------------------------------------------------------]
[--------c-------------------------------------------------------]
[------------d---------------------------------------------------]
etc...

*/

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE typename spectral<T, SpectrumSize, Overlap>::sample_t 
spectral<T, SpectrumSize, Overlap>::generate()
{
  if (overlap_count_ == block_size)
  {
    fill_spectrum();
    overlap_add();
    overlap_count_ = 0;
  }

  sample_t& read = output_[read_idx_];
  sample_t out   = read;
  read = 0;
  read_idx_ = (read_idx_+1)&(SpectrumSize-1);
  ++overlap_count_;
  
  return out;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::generate(array<T> out)
{
  VASSERT(out.size() == block_size, "Output array for spectral::generate(array) must be exactly block_size long");

  fill_spectrum();
  overlap_add();

  array<T> out_block(output_.data() + read_idx_, block_size);
  out_block.copy_to(out);
  out_block.fill(0);
  read_idx_ = (read_idx_+block_size)&(SpectrumSize-1);
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE size_t spectral<T, SpectrumSize, Overlap>::get_read_head(size_t idx) const
{
  return read_idx_;
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::fill_spectrum()
{
  // gained a better understanding of glitching.
  // it stems from generating complex numbers for spectrum_
  // from band magnitudes that are "too small,"
  // which I think is relative to SpectrumSize.
  // magnitude needs to be scaled up for a band quite a bit when setting it externally
  // in order to generate a time-domain signal that is at the level I expect.
  // experiment with where this amplification occurs.
  // is it better applied here, or after generation of a time-domain frame?
  // users of this class should be able to set a [0,1] magnitude in a band
  // and get comparably leveled audio at the output.
  static constexpr T mag_zero  = cast<T>(0);
  static constexpr T mag_min   = cast<T>(1.f/SpectrumSize);
  static constexpr T mag_scale = cast<T>(static_cast<analog_t>(SpectrumSize)/32.f);
  // DC component
  spectrum_[0].set_complex(0,0);
  const size_t max_band = spectrum_.size() - 1;
  for (size_t i = 1; i < max_band; ++i)
  {
    const frequency_band& band = bands_[i];
    T m = band.magnitude() > mag_min ? band.magnitude() * mag_scale : mag_zero;
    T s = 1 - 2*((i&1)&phase_flip_);
    complex_t cmplx = band.to_complex();
    cmplx.scale(s*m);
    // //if (m > mag_zero)
    // {
    //   if constexpr (FlipOddPhases)
    //   {
    //     T s = i&1 ? -1 : 1;
    //     cmplx = band.to_complex();
    //     cmplx.scale(s*m);
    //   }
    //   else
    //   {
    //     cmplx = band.to_complex();
    //     cmplx.scale(m);
    //   }
    // }
    spectrum_[i] = cmplx;
  }
  phase_flip_ = 1 - phase_flip_;
}

template <typename T, size_t SpectrumSize, size_t Overlap> 
VESSL_INLINE void spectral<T, SpectrumSize, Overlap>::overlap_add()
{
  fft_.inverse(spectrum_, signal_);
  signal_.multiply(window_);

  const size_t head_size = output_.size() - read_idx_;
  array<T> sig_head(signal_.data(), head_size);
  array<T> out_head(output_.data() + read_idx_, head_size);
  out_head.add(sig_head);

  // copy tail only if we have one
  if (const size_t tail_size = read_idx_)
  {
    array<T> sig_tail(signal_.data() + head_size, tail_size);
    array<T> out_tail(output_.data(), tail_size);  
    out_tail.add(sig_tail);
  }
}
} // namespace generators
} // namespace vessl