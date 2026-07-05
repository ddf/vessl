#pragma once

namespace vessl
{
namespace processors
{
template<typename T>
VESSL_INLINE T slew<T>::process(const T& v)
{
  binary_t is_rise = v > params_.output.value+eps_;
  binary_t is_fall = v < params_.output.value-eps_;
  if (is_rise)
  {
    params_.output.value = math::min(v, params_.output.value + params_.rise.value*dt_);
  }
  else if (is_fall)
  {
    params_.output.value = math::max(v, params_.output.value - params_.fall.value*dt_);
  }
  else
  {
    params_.output.value = v;
  }
  params_.rising.value = is_rise;
  params_.falling.value = is_fall;
  return params_.output.value;
}

template<typename T>
VESSL_INLINE T delay<T>::process(const T& in)
{
  delay_in_samples_ = params_.time.value.samples;
  // delay time in samples
  analog_t dts = math::constrain<analog_t>(delay_in_samples_, 0.f, cast<analog_t>(buffer_.size()-1));
  analog_t s = buffer_.readf(dts);
  analog_t fbk = math::constrain<analog_t>(cast<analog_t>(feedback()), -1.0, 1.0);
  buffer_.write(in + s * fbk);
  return s;
}

template<typename T>
template<time::mode TimeMode>
VESSL_INLINE void delay<T>::process(array<T> input, array<T> output)
{
  if (TimeMode == time::mode::snap)
  {
    processor<T>::process(input, output);
  }

  if (TimeMode == time::mode::slew)
  {
    auto r = input.make_reader();
    auto w = output.make_writer();
    analog_t dst = dt_ * 10.0f;
    while (r && w)
    {
      T in = r.read();
      delay_in_samples_ = math::lerp(delay_in_samples_, params_.time.value.samples, dst);
      // delay time in samples
      analog_t dts = math::constrain<analog_t>(delay_in_samples_, 0.f, cast<analog_t>(buffer_.size()-1));
      analog_t wet = buffer_.readf(dts);
      analog_t fbk = math::constrain<analog_t>(cast<analog_t>(feedback()), -1.0, 1.0);
      buffer_.write(in + wet*fbk);
      w << wet;
    }
  }
  
  if (TimeMode == time::mode::fade)
  {
    auto r = input.make_reader();
    auto w = output.make_writer();
    
    analog_t fade = 0;
    analog_t fadeInc = 1.0f / input.size();
    // smooth time parameter to prevent crunchiness when it is noisy or changes by large amounts
    analog_t targetSampleDelay = params_.time.value.samples;
    // delay time in samples
    analog_t sza = cast<analog_t>(buffer_.size()-1);
    analog_t fts = math::constrain<analog_t>(delay_in_samples_, 0.0, sza);
    analog_t tts = math::constrain<analog_t>(targetSampleDelay, 0.0, sza);
    analog_t fbk = math::constrain<analog_t>(cast<analog_t>(feedback()), -1.0, 1.0);
    
    while (r && w)
    {
      T in = r.read();
      T wet = (1.0f - fade) * buffer_.readf(fts) + fade * buffer_.readf(tts);
      buffer_.write(in + wet*fbk);
      w << wet;
      fade += fadeInc;
    }
      
    delay_in_samples_ = targetSampleDelay;
  }
}

template<typename T>
VESSL_INLINE T follow<T>::process(const T& in) 
{
  writer_.write(in);
  if (writer_.is_full())
  {
    previous_ = current_;
    current_ = T(0);
    auto r = window_.make_reader();
    while (r)
    {
      current_ *= delta_;
      current_ += (1.0 - delta_)*math::abs(r.read());
    }
    writer_ = window_.make_writer();
  }

  analog_t t = 1.0 - cast<analog_t>(writer_.available()) / cast<analog_t>(window_.size());
  return previous_ + (current_ - previous_) * t;
}

template <typename T>
typename processor<T>::output_t peak_meter<T>::process(const typename processor<T>::input_t &in)
{  
  T peak = math::abs(in);
  T error = peak - params_.peak.value;
  params_.peak.value += (error > T(0)) ? T(0.05) : T(0.00002) * error;
  return params_.peak.value;
}

template<typename T>
VESSL_INLINE T limiter<T>::process(const T& in)
{
  T pre = in*params_.pre_gain.value.to_scale();
  T peak = peak_meter<T>::process(in);
  analog_t gain = (peak <= 1.0 ? 1.0 : 1.0 / peak);
  // DaisySP returns this, which sounds better for how I typically use this.
  return sample::softlimit(pre*gain*0.7);
  // stmlib returns this, which clips more easily, but is faster.
  //return pre*gain*0.8;
}

template<typename T>
VESSL_INLINE T freeze<T>::generate() 
{
  freeze_delay_ = cast<analog_t>(position());
  freeze_size_  = params_.duration.value.samples;
  analog_t sampleDelay = freeze_delay_ + (1.0-phase_)*freeze_size_;
  phase_ = math::wrap01(phase_ + rate() / freeze_size_);
  return delay_line_.readf(sampleDelay);
}

template<typename T>
VESSL_INLINE T freeze<T>::process(const T& in) 
{
  binary_t is_enabled = cast<binary_t>(enabled());
  crossfade_ = (is_enabled ? 1.0 : 0.0);
  analog_t wet_level = crossfade_.value;
  T wet = generate();
  if (!is_enabled)
  {
    delay_line_.write(in);
  }
  return wet_level*wet + (1.0 - wet_level)*in;
}

template<typename T>
template<time::mode TimeMode, bool UseInput>
VESSL_INLINE void freeze<T>::proc_gen(array<T> input, array<T> output) 
{
  typename array<T>::reader r(input);
  typename array<T>::writer w(output);

  if (TimeMode == time::mode::snap)
  {
    if (UseInput)
    {
      while (w)
      {
        w << process(r.read());
      }
    }
    else
    {
      while (w)
      {
        w << generate();
      }
    }
  }

  if (TimeMode == time::mode::slew)
  {
    analog_t fd = params_.position.value;
    analog_t fs = params_.duration.value.samples;
    analog_t rt = params_.rate.value;
    analog_t st = dt_;
    while (w)
    {
      freeze_delay_ = math::lerp(freeze_delay_, fd, st*20);
      freeze_size_ = math::lerp(freeze_size_, fs, st*20);
      read_rate_ = math::lerp(read_rate_, rt, st*10);
      analog_t sample_delay = freeze_delay_ + (1.0 - phase_)*freeze_size_;
      phase_ = math::wrap01(phase_ + read_rate_/freeze_size_);
      T wet = delay_line_.readf(sample_delay);

      if (UseInput)
      {
        binary_t is_enabled = params_.enabled.value;
        crossfade_ = (is_enabled ? 1.0 : 0.0);
        analog_t wet_level = crossfade_.value;
        T in = r.read();
        if (!is_enabled)
        {
          delay_line_.write(in);
        }
        w << wet*wet_level + in*(1.0 - wet_level);
      }
      else
      {
        w << wet;
      }
    }
  }

  if (TimeMode == time::mode::fade)
  {
    analog_t fd0 = freeze_delay_,   fd1 = params_.position.value;
    analog_t fs0 = freeze_size_,    fs1 = params_.duration.value.samples;
    analog_t fade = 0, fadeInc = 1.0 / output.size();
    analog_t r0 = read_rate_, r1 = params_.rate.value;
    analog_t p0 = phase_, dp0 = r0/fs0, dp1 = r1/fs1;
    while (w)
    {
      analog_t sd0 = fd0 + fs0*(1.0f-p0);
      analog_t sd1 = fd1 + fs1*(1.0f-phase_);
      T wet = (1.0 - fade)*delay_line_.readf(sd0) + fade*delay_line_.readf(sd1);
      phase_ = math::wrap01(phase_ + dp1);
      p0 = math::wrap01(p0 + dp0);
      fade += fadeInc;

      if (UseInput)
      {
        binary_t is_enabled = params_.enabled.value;
        crossfade_ = (is_enabled ? 1.0f : 0.0f);
        analog_t wet_level = crossfade_.value;
        T in = r.read();
        if (!is_enabled)
        {
          delay_line_.write(in);
        }
        w << wet_level*wet + (1.0 - wet_level)*in;
      }
      else
      {
        w << wet;
      }
    }
    
    freeze_delay_ = fd1;
    freeze_size_ = fs1;
    read_rate_ = r1;
  }
}

template<typename T, uint32_t MaxBits>
VESSL_INLINE T bitcrush<T, MaxBits>::process(const T& in)
{
  rate_alpha_ += math::max<analog_t>(1.0, cast<analog_t>(rate()))*dt_;
  if (rate_alpha_ >= 1)
  {
    rate_alpha_ -= 1;
    curr_sample_ = math::lerp(prev_input_, in, rate_alpha_);
  }

  analog_t bd = math::constrain<analog_t>(cast<analog_t>(depth()), 2.0, MaxBits);
  analog_t scalar = math::pow<analog_t>(2.0, bd) - 1;
  T val = math::round(curr_sample_*scalar);
  if (mangle())
  {
    val = math::xore(val, math::round(prev_input_*scalar));
  }
  prev_input_ = in;
  return val * (1.0 / scalar);
}
} // namespace processors
} // namespace vessl