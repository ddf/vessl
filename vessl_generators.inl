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

template <typename T, size_t SpectrumSize, size_t Overlap>
spectral<T, SpectrumSize, Overlap>::spectral(data &data, analog_t sample_rate)
  : unit_generator<T>()
  , fft_(SpectrumSize)
  , frequencies_(data.frequencies.data(), data.frequencies.size())
  , spectrum_(data.spectrum.data(), data.spectrum.size())
  , signal_(data.signal.data(), data.signal.size())
  , window_(data.window.data(), data.window.size())
  , buffer_(data.buffer.data(), data.buffer.size())
  , read_idx_(0)
  , gen_idx_(0)
  , gen_inc_(SpectrumSize/Overlap)
  , phase_shift_(0)
{
  constexpr size_t bands = SpectrumSize/2;
  VASSERT(data.spectrum.frequencies.size() == bands, "Invalid frequency bands size");
  VASSERT(data.spectrum.size() >= bands, "Invalid spectrum size");
  VASSERT(data.window.size() == SpectrumSize, "Invalid window size");
  VASSERT(data.signal.size() == SpectrumSize, "Invalid signal size");
  VASSERT(data.buffer.size() >= bands, "Invalid buffer size");
  spectrum_[0] = { 0, 0 };
  for (int i = 0; i < bands; ++i)
  {
    frequencies_[i].magnitude = 0;
    frequencies_[i].phase = math::random::u32();
  }
}

template <typename T, size_t SpectrumSize, size_t Overlap>
VESSL_INLINE typename spectral<T, SpectrumSize, Overlap>::sample_t spectral<T, SpectrumSize, Overlap>::generate()
{
  if (gen_idx_ == read_idx_)
  {
    frequency_band* bands = frequencies_.data();
    spectrum_.fill(complex_t(0.));
    complex_t* spectrum = spectrum_.data();
    for (int i = 0; i < SpectrumSize/2 - 1; ++i)
    {
      size_t si = i+1;
      phase_t ps = 0; // phase_shift_*(si%2);
      spectrum[si].set_polar(bands[i].magnitude, bands[i].phase + ps);
    }
    fft_.inverse(spectrum_, signal_);
    //data_.signal.multiply(data_.window);
    for (int i = 0; i < SpectrumSize; ++i)
    {
      buffer_.overdub(signal_[i]*window_[i]);
      //data_.buffer.write(r.read());
    }
    size_t widx = buffer_.get_write_index();
    gen_idx_ = widx > gen_inc_ ? widx - gen_idx_ - 1 : buffer_.size() - 1 - (gen_inc_ - widx);
    phase_shift_ = (phase_shift_ == 0 ? phase_180 : 0);
  }

  sample_t sample = buffer_.data()[read_idx_];
  buffer_.data()[read_idx_] = 0.f;
  read_idx_ = (read_idx_ + 1) % buffer_.size();
  return sample;
}
} // namespace generators
} // namespace vessl