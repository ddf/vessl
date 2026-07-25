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

template <typename T, size_t SpectrumSize>
spectral<T, SpectrumSize>::spectral(data &data, analog_t sample_rate)
  : unit_generator<T>()
  , fft_(SpectrumSize)
  , frequencies_(data.frequencies.data(), data.frequencies.size())
  , spectrum_(data.spectrum.data(), data.spectrum.size())
  , signal_a_(data.signal.data(), SpectrumSize)
  , signal_b_(data.signal.data() + SpectrumSize, SpectrumSize)
  , window_(data.window.data(), data.window.size())
  , read_idx_a_(SpectrumSize)
  , read_idx_b_(SpectrumSize/2)
{
  constexpr size_t bands = SpectrumSize/2;
  VASSERT(data.frequencies.size() == bands, "Invalid frequency bands size");
  VASSERT(data.spectrum.size() == bands, "Invalid spectrum size");
  VASSERT(data.window.size() == SpectrumSize, "Invalid window size");
  VASSERT(data.signal.size() == SpectrumSize*Overlap, "Invalid signal size");

  for (int i = 0; i < bands; ++i)
  {
    frequency_band& band = frequencies_[i];
    band.magnitude = 0;
    band.phase = math::random::u32();
  }
  signal_a_.fill(0);
  signal_b_.fill(0);
}

template <typename T, size_t SpectrumSize>
VESSL_INLINE typename spectral<T, SpectrumSize>::sample_t spectral<T, SpectrumSize>::generate()
{
  if (read_idx_a_ == SpectrumSize)
  {
    fill_spectrum<false>();
    fft_.inverse(spectrum_, signal_a_);
    read_idx_a_ = 0;
  }
  
  if (read_idx_b_ == SpectrumSize)
  {
    fill_spectrum<true>();
    fft_.inverse(spectrum_, signal_b_);
    read_idx_b_ = 0;
  }
  
  sample_t out = 0;
  out += signal_a_[read_idx_a_] * window_[read_idx_a_];
  ++read_idx_a_;
  out += signal_b_[read_idx_b_] * window_[read_idx_b_];
  ++read_idx_b_;
  
  return out;
}

template <typename T, size_t SpectrumSize>
VESSL_INLINE size_t spectral<T, SpectrumSize>::get_read_head(size_t idx) const
{
  return idx == 0 ? read_idx_a_ : read_idx_b_;
}

template <typename T, size_t SpectrumSize>
template<bool ShiftOddPhases>
VESSL_INLINE void spectral<T, SpectrumSize>::fill_spectrum()
{
  spectrum_[0].set_complex(0,0);
  for (int i = 1; i < spectrum_.size(); ++i)
  {
    frequency_band& band = frequencies_[i-1];
    T m = band.magnitude;
    // for reasons I do not understand, this only works correctly if z is fully assigned.
    // initializing z to band.phase and then incrementing it by phase_180 
    // when ShiftOddPhases is true creates artifacts.
    // similarly, if band.magnitude is zero and we don't use phase_zero for the phase,
    // we also generate artifacts.
    phase_t z;
    if constexpr (ShiftOddPhases)
    {
      z = m>0 ? (i&1 ? band.phase + phase_180 : band.phase) : phase_zero; 
    }
    else
    {
      z = m>0 ? band.phase : phase_zero;
    }
    spectrum_[i].set_polar(m, z);
  }
}
} // namespace generators
} // namespace vessl