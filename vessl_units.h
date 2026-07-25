#pragma once

#include "vessl.h"
#include <algorithm> // for swap

namespace vessl
{
namespace generators
{
// generates stepped analog noise at rate values per second.
template<typename T, typename N>
class noise : public unit_generator<T>, protected plist<1>
{
public:
  N noise_source;
    
  explicit noise(analog_t sample_rate = 1)
  : unit_generator<T>(), noise_source(sample_rate), dt_(1.0f/sample_rate), step_(0)
  {
    value_ = noise_source();
    next_ = noise_source();
    params_.rate.value = sample_rate;
  }
    
  void set_sample_rate(float sample_rate) override { dt_ = 1.0f / sample_rate;}
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }

  [[nodiscard]] parameter rate() const { return params_.rate("rate", 'r'); }

  // generates stepped noise in the range [0,1] at the given rate
  T generate() override;

  // smooths the stepped noise with the given easing
  template<typename E>
  T generate();

protected:
  [[nodiscard]] parameter element_at(size_t index) const override { parameter p[num] = { rate() }; return p[index]; }

private:
  struct
  {
    analog_p rate;
  } params_;
  analog_t dt_;
  analog_t value_;
  analog_t next_;
  analog_t step_;
};

// unit that generates a linear ramp from one value to another over a duration of seconds
// @todo implement easings above and add that as a template parameter
template<typename T>
class ramp : public unit_generator<T>, protected plist<4>
{
public:
  explicit ramp(analog_t sample_rate, analog_t duration_in_seconds = 0, T from_value = T(0), T to_value = T(0))
  : unit_generator<T>(), dt_(1.0f/sample_rate), t_(0)
  {
    params_.from.value = from_value;
    params_.to.value = to_value;
    params_.duration.value = duration_in_seconds;
    params_.eor.value = false;
  }
    
  void set_sample_rate(float sample_rate) override { dt_ = 1.0f / sample_rate;}
  const parameter_list& parameters() const override { return *this; }

  // ins
  parameter from() const { return params_.from("from", 'f'); }
  parameter to() const { return params_.to("to", 't'); }
  parameter duration() const { return params_.duration("duration", 'd'); }

  // outs
  parameter eor() const { return params_.eor("eor", 'e'); }
  // could also add t as an out.

  binary_t is_active() const { return !params_.eor.value; }
  void trigger();
  T generate() override;
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = { from(), to(), duration(), eor() }; return p[index];
  }
    
private:
  struct
  {
    param<T> from, to;
    analog_p duration;
    binary_p eor;
  } params_;
  analog_t dt_;
  analog_t t_;
};

// generates an envelope that begins and ends at zero, with some number of stages leading up to a final stage.
// each stage of an envelope is defined by a target value, a duration to reach it,
// and whether the envelope should hold the stage value until it is triggered again.
// @todo figure out how the heck to return a full description and parameter lists for envelopes.
template<typename T>
class envelope : public unit_generator<T>, protected plist<2>
{
public:
  // @todo move this functionality into ramp and make envelope a series of ramps
  class stage final : public unit_generator<T>, protected plist<5>
  {
  public:
    explicit stage(analog_t sample_rate) : unit_generator<T>(), begin_(0), dt_(1.0f/sample_rate) { reset(); }
      
    void set_sample_rate(float sample_rate) override { dt_ = 1.0f / sample_rate;}
    const parameter_list& parameters() const override { return *this; }

    parameter target() const { return params_.target("target", 't'); }
    parameter duration() const { return params_.duration("duration", 'd'); }
    parameter active() const { return params_.active("active", 'a'); }
    parameter eos() const { return params_.eos("eos", 'e'); }
    // current value of the stage
    parameter value() const { return params_.output("value", 'v'); }

    void start(T from_value) { begin_ = from_value; params_.output.value = from_value; time_ = -dt_; params_.active.value = true; params_.eos.value = false; }
    void reset() { params_.active.value = false; params_.eos.value = false; time_ = -dt_; params_.output.value = 0; }
      
    template<typename E>
    T generate() { return active() ? step<E>() : params_.output.value; }
    T generate() override { return generate<math::easing::linear>(); }
      
  protected:
    parameter element_at(size_t index) const override
    {
      parameter p[num] = { target(), duration(), active(), eos(), value() };
      return p[index];
    }
      
  private:
    template<typename E>
    T step();
      
    struct
    {
      param<T> target;
      analog_p duration;
      binary_p active;
      binary_p eos;
      param<T> output;
    } params_;
    T begin_; // value the stage started with
    // where we are in the stage
    analog_t time_;
    analog_t dt_;
  };
    
  void set_sample_rate(float sample_rate) override;
    
  stage& get_stage(size_t idx) { return idx == stages_.size() ? final_ : stages_[idx]; }
  const stage& get_stage(size_t idx) const { return idx == stages_.size() ? final_ : stages_[idx]; }
  size_t get_stage_count() const { return stages_.size() + 1; }

  stage& current_stage() { return get_stage(stage_idx_); }
  const stage& current_stage() const { return get_stage(stage_idx_);}
  stage& final_stage() { return final_; }
  const stage& final_stage() const { return final_; }
    
  const parameter_list& parameters() const override { return *this; }

  parameter value() const { return current_stage().value(); }
  parameter eoc() const { return params_.eoc("eoc", 'e'); }

  // make this a parameter we check in generate?
  virtual void trigger();

  template<typename E>
  T generate();
  T generate() override { return generate<math::easing::linear>(); }

protected:
  envelope(stage* stage_data, size_t stage_data_size, analog_t sample_rate) : unit_generator<T>()
  , stages_(stage_data, stage_data_size), stage_idx_(0), final_(sample_rate) {}
    
  parameter element_at(size_t index) const override { parameter p[num] = { value(), eoc() }; return p[index]; }
  void start_stage(size_t idx, T from_value) { get_stage(idx).start(from_value); stage_idx_ = idx; }
    
  // by default, stages advance automatically when their eos goes high.
  // subclasses can override this behavior per stage
  // to enable advancing to the next stage before it is finished,
  // or holding a stage for some period of time.
  virtual binary_t should_advance(size_t current_stage_idx) { return get_stage(current_stage_idx).eos().template read<bool>(); }
  
private:
  struct
  {
    binary_p eoc;
  } params_;
  array<stage> stages_;
  size_t stage_idx_;
  stage final_;
};

template<typename T>
class ad : public envelope<T>
{
public:
  ad(analog_t attack_duration, analog_t decay_duration, analog_t sample_rate) 
  : envelope<T>(&attack_stage_, 1, sample_rate), attack_stage_(sample_rate)
  {
    attack().target() = T(1); attack().duration() = attack_duration; decay().duration() = decay_duration;
  }
    
  typename envelope<T>::stage& attack() { return attack_stage_; }
  typename envelope<T>::stage& decay() { return envelope<T>::final_stage(); }
  using envelope<T>::eoc;

  using envelope<T>::trigger;
  using envelope<T>::generate;
    
protected:
  using envelope<T>::element_at;
    
private:
  typename envelope<T>::stage attack_stage_;
};

template<typename T>
class asr : public ad<T>
{
public:
  asr(analog_t attack_duration, analog_t decay_duration, analog_t sample_rate, T trigger_threshold = T(0))
  : ad<T>(attack_duration, decay_duration, sample_rate)
  , trig_threshold_(trigger_threshold), gate_on_(false) 
  {}

  using ad<T>::attack;
  typename envelope<T>::stage& release() { return envelope<T>::final_stage(); }
  using ad<T>::eoc; 

  void gate(T value);
  void gate(binary_t on) { gate(on ? T(1) : T(0)); }
  void trigger() override { attack().target() = T(1); ad<T>::trigger(); gate_on_ = false; }
  using envelope<T>::generate;

protected:
  binary_t should_advance(size_t current_stage_idx) override
  {
    return ad<T>::should_advance(current_stage_idx) && !gate_on_;
  }
    
private:
  T trig_threshold_;
  binary_t gate_on_;
};

template<typename T>
class adsr : public envelope<T>
{
public:
  adsr(analog_t attack_duration, analog_t decay_duration, analog_t sustain_level, analog_t release_duration
    , analog_t sample_rate)
  : envelope<T>(&attack_stage_, 2, sample_rate), attack_stage_(sample_rate), decay_stage_(sample_rate), gate_on_(false)
  {
    attack_stage_.duration() = attack_duration;
    attack_stage_.target() = T(1);
    decay_stage_.duration() = decay_duration;
    decay_stage_.target() = sustain_level;
    envelope<T>::final_stage().duration() = release_duration;
  }

  typename envelope<T>::stage& attack() { return attack_stage_; }
  typename envelope<T>::stage& decay() { return decay_stage_; }
  parameter& sustain() { return decay().target(); }
  typename envelope<T>::stage& release() { return envelope<T>::finalStage(); }
  using envelope<T>::eoc;

  // should we jump to the release stage if the gate goes off before we start sustaining??
  void gate(binary_t on)
  {
    if (on && !gate_on_)
    {
      envelope<T>::trigger();
    }
    else if (gate_on_ && !on)
    {
      envelope<T>::start_stage(2, envelope<T>::current_stage().value().template read<T>());
    }
    gate_on_ = on;
  }
  void trigger() override { gate_on_ = false; envelope<T>::trigger(); }
  using envelope<T>::generate;

protected:
  binary_t should_advance(size_t current_stage_idx) override
  {
    return current_stage_idx == 1 ? decay_stage_.eos() && !gate_on_ : envelope<T>::should_advance(current_stage_idx);
  }
    
private:
  typename envelope<T>::stage attack_stage_;
  typename envelope<T>::stage decay_stage_;
  binary_t gate_on_;
};

// note: W must implement waveform<T>
template<class W>
class oscil final : public unit_generator<typename W::sample_t>, protected plist<4>
{
public:
  W waveform;
  using sample_t = typename W::sample_t;
    
  oscil() : unit_generator<sample_t>(), phase_(phase_zero), dt_(phase_zero) { params_.fhz.value = 440.0; }

  template<typename... Ts>
  explicit oscil(analog_t sample_rate, analog_t freq_in_hz, Ts... wargs)
  : unit_generator<sample_t>(), waveform(wargs...), phase_(phase_zero), dt_(cast<phase_t>(1.0f/sample_rate))
  {
    params_.fhz.value = freq_in_hz;
  }
    
  void set_sample_rate(float sample_rate) override { dt_ = cast<phase_t>(1.0f/sample_rate); }
    
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
    
  // frequency in Hz without FM applied
  [[nodiscard]] parameter fhz() const { return params_.fhz("frequency", 'f'); }
  // linear frequency modulation
  [[nodiscard]] parameter fm_lin() const { return params_.fm_lin("fm (lin)", 'l'); }
  // v/oct (exponential) frequency modulation
  [[nodiscard]] parameter fm_exp() const { return params_.fm_exp("fm (v/oct)", 'v'); }
  // phase modulation
  [[nodiscard]] parameter pm() const { return params_.pm("phase mod", 'p'); }

  sample_t generate() override;
  void generate(sink<sample_t>& dest);

  [[nodiscard]] phase_t phase() const { return phase_; }
  [[nodiscard]] phase_t inc() const { return dt_; }
  void reset() { phase_ = phase_zero; }
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] { fhz(), fm_lin(), fm_exp(), pm() };
    return p[index];
  }
    
private:
  struct
  {
    analog_p fhz;
    analog_p fm_lin;
    analog_p fm_exp;
    phase_p  pm;
  } params_;
  phase_t  phase_;
  phase_t  dt_;
};

template<typename T>
class clock final : public unit_generator<T>
  , protected plist<1>
  , time::clockable
{
public:
  clock(analog_t sample_rate, period_t sample_period_min, period_t sample_period_max, analog_t bpm = 60);

  [[nodiscard]] const parameter_list & parameters() const override { return *this; }
  
  VESSL_INLINE void tap() { clockable::clock(); }
  VESSL_INLINE void tap(period_t sample_delay) { clockable::clock(sample_delay); }
  
  using clockable::is_clocked;
  
  // tempo of the clock expressed as a duration
  // so that it can be converted to BPM, sample period, or frequency as desired.
  [[nodiscard]] VESSL_INLINE parameter tempo() const;
  
  T generate() override;

protected:
  [[nodiscard]] parameter element_at(size_t index) const override;
  
private:
  phase_t dt_;
  phase_t phase_;
  sample::waves::unipolar::square<T> pulse_;
};

// synthesize audio from a frequency-domain spectrum.
// SpectrumSize is the size of the FFT to use.
// This will determine the length of the blocks of time-domain signal
// that are overlapped to produce this generator's output.
// The number of frequency bands available will be SpectrumSize/2 - 1.
template<typename T, size_t SpectrumSize, size_t Overlap = 2>
class spectral : public unit_generator<T>, plist<0>
{
public:
  using fft_t = transform::fft<T>;
  using sample_t = T;
  using complex_t = transform::complex<T>;
  
  struct frequency_band
  {
    T magnitude = 0;
    phase_t phase = 0;
  };
  
  // data.frequencies must have length equal to SpectrumSize/2
  // data.spectrum must have a length equal to SpectrumSize/2
  // data.signal must have a length equal to SpectrumSize
  // data.window must have a length equal to SpectrumSize
  // data.buffer must have a length greater than or equal to SpectrumSize/2
  struct data
  {
    array<frequency_band> frequencies;
    array<complex_t> spectrum;
    array<sample_t> signal;
    array<sample_t> window;
    sample::ring_buffer<sample_t> buffer;
  };
  
  spectral(data& data, analog_t sample_rate);
  
  [[nodiscard]] const parameter_list & parameters() const override { return *this; }
  [[nodiscard]] VESSL_INLINE frequency_band& get_band(size_t index) { return frequencies_[index]; }
  sample_t generate() override;

protected:
  fft_t fft_;
  array<frequency_band> frequencies_;
  array<complex_t> spectrum_;
  array<sample_t> signal_;
  array<sample_t> window_;
  sample::ring_buffer<sample_t> buffer_;
  size_t read_idx_;
  size_t gen_idx_;
  size_t gen_inc_;
  phase_t phase_shift_;
};

} // namespace generators
  
namespace processors
{
template<typename T>
class slew : public unit_processor<T>, protected plist<5>
{
public:
  // note: choice of epsilon will depend on the amount of noise in the signal to be slewed.
  // the default value was chosen based on testing with an OWL module's audio input.
  slew(analog_t sample_rate, analog_t rise_rate, analog_t fall_rate, T initial_value = T(0), T epsilon = math::epsilon<T>()*1000)
  : unit_processor<T>(), eps_(epsilon), dt_(1.0f/sample_rate)
  {
    params_.rise.value = rise_rate; params_.fall.value = fall_rate; params_.output.value = initial_value;
  }
    
  void set_sample_rate(float sample_rate) override { dt_ = 1.0f / sample_rate; }
  const parameter_list& parameters() const override { return *this; }

  parameter rise() const { return params_.rise("rise", 'a'); }
  parameter fall() const { return params_.fall("fall", 'd'); }
  parameter rising() const { return params_.rising("rising", 'r'); }
  parameter falling() const { return params_.falling("falling", 'f'); }
  parameter value() const { return params_.output("value", 'v'); }
    
  T process(const T& v) override;
  using processor<T>::process;
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = { rise(), fall(), rising(), falling(), value() };
    return p[index];
  }
    
private:
  struct
  {
    analog_p rise;
    analog_p fall;
    binary_p rising;
    binary_p falling;
    param<T> output;
  } params_;
  T eps_;
  analog_t dt_;
};
  
template<typename T>
class delay : public unit_processor<T>, protected plist<2>
{
public:
  delay(array<T> delay_buffer, analog_t sample_rate, analog_t delay_in_seconds = 0, analog_t feedback_amount = 0)
  : unit_processor<T>(), buffer_(delay_buffer.data(), delay_buffer.size()), dt_(1.0f/sample_rate)
  {
    params_.time.value = duration_t::from_seconds(delay_in_seconds, sample_rate);
    delay_in_samples_ = params_.time.value.samples;
    params_.feedback.value = feedback_amount;
  }
    
  void set_sample_rate(float sampleRate) override { dt_ = 1.0f / sampleRate; }
  const parameter_list& parameters() const override { return *this; }

  sample::delay_line<T>& buffer() { return buffer_; }
  const sample::delay_line<T>& buffer() const { return buffer_; }

  /// delay time expressed as vessl::duration (i.e. samples), can be set using an analog_t
  parameter time() const { return params_.time("time", 't'); }
  /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
  parameter feedback() const { return params_.feedback("feedback", 'f'); }
    
  T process(const T& in) override;

  void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
  template<time::mode TimeMode = time::mode::slew>
  void process(array<T> input, array<T> output);
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = { time(), feedback() }; 
    return p[index];
  }
    
  sample::delay_line<T> buffer_;
  analog_t delay_in_samples_;
  analog_t dt_;

private:
  struct
  {
    duration_p time;
    analog_p feedback;
  } params_;
};

template<typename T>
class follow : public unit_processor<T>, protected plist<1>
{
public:
  follow(array<T> window_array, analog_t sample_rate, analog_t response_time_in_seconds)
  : unit_processor<T>(), writer_(window_array), window_(window_array)
  , delta_(math::exp(-1.0 / (sample_rate*response_time_in_seconds)))
  , previous_(0), current_(0)
  {
    params_.response.value = response_time_in_seconds;
  }
    
  void set_sample_rate(float sample_rate) override { delta_ = math::exp(-1.0 / (sample_rate*params_.response.value)); }
  const parameter_list& parameters() const override { return *this; }
    
  parameter response() const { return params_.response("response time", 'r'); }

  T process(const T& in) override;

  using processor<T>::process;
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = { response() }; return p[index];
  }
    
private:
  struct
  {
    // @todo actually use this parameter
    analog_p response;
  } params_;
    
  typename
  array<T>::writer writer_;
  array<T> window_;
  analog_t delta_;
  T previous_;
  T current_;
};

template<typename T>
class peak_meter : public unit_processor<T>, protected plist<1>
{
public:
  
  [[nodiscard]] const parameter_list & parameters() const override { return *this; }
  [[nodiscard]] parameter peak() const { return params_.peak("peak", 'k'); }
  typename processor<T>::output_t process(const typename processor<T>::input_t &in) override;
  
protected:
  [[nodiscard]] VESSL_INLINE parameter element_at(size_t index) const override
  {
    return index == 0 ? peak() : parameter::none();
  }
  
private:
  struct 
  {
    param<T> peak;
  } params_;
};

// A simple peak limiter adapted from pinchenettes/stmlib via DaisySP
template<typename T>
class limiter : public peak_meter<T>
{
public:
  explicit limiter(gain_t pre_gain = gain_t::from_decibels(0))
  {
    params_.pre_gain.value = pre_gain;
  }

  [[nodiscard]] const parameter_list& parameters() const override { return *this; }

  [[nodiscard]] parameter pre_gain() const { return params_.pre_gain("pre-gain", 'g'); }
  using peak_meter<T>::peak;

  T process(const T& in) override;
  using unit_processor<T>::process;
    
protected:
  [[nodiscard]] VESSL_INLINE size_t size() const override { return 2; }
  [[nodiscard]] parameter element_at(size_t index) const override
  {
    switch (index)
    {
    case 0: return pre_gain();
    case 1: return peak();
    default: return parameter::none();
    }
  }
    
private:
  struct
  {
    gain_p pre_gain;
  } params_;
};
  
// when used as a processor, will write incoming audio to the delayline
// and output the incoming signal if freeze is not engaged.
// when freeze is engaged, it will ignore the incoming signal
// and generate audio from previously recorded input.
//
// when used as a generator, it will always generate using the
// contents of the delayline, enabling it to be used to "freeze"
// audio that has been recorded elsewhere (e.g. by a delay).
template<typename T>
class freeze : public unit, public processor<T>, public generator<T>, protected plist<4>
{
public:
  explicit freeze(array<T> freeze_buffer, analog_t sample_rate) : unit()
  , delay_line_(freeze_buffer.data(), freeze_buffer.size())
  , phase_(0), crossfade_(0.75f)
  , freeze_delay_(0), freeze_size_(0)
  , read_rate_(1), dt_(1.0f/sample_rate)
  {
    params_.duration.value.samples = delay_line_.size()-1;
    freeze_size_ = params_.duration.value.samples;
    params_.rate.value = 1.0;
  }
    
  void set_sample_rate(analog_t sr) override { dt_ = 1.0f/sr; }
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
    
  sample::delay_line<T>& buffer() { return delay_line_; }
  const sample::delay_line<T>& buffer() const { return delay_line_; }
    
  parameter enabled() const { return params_.enabled("enabled", 'e'); }
  // end of the freeze loop in samples relative to the most recently recorded sample
  parameter position() const { return params_.position("position", 'p'); }
  // size of the freeze loop as a duration (samples).
  // the beginning of the freeze loop, when played forward, will be position + size.
  parameter duration() const { return params_.duration("duration", 'd'); }
  // rate of playback when enabled, can be negative to play in reverse
  parameter rate() const { return params_.rate("rate", 'r'); }
  // should this be a parameter? there's not much gained by it.
  analog_t phase() const { return phase_; }
  // reset the phase to zero (argument for making it a parameter?)
  void reset() { phase_ = 0; }

  T generate() override;

  T process(const T& in) override;

  VESSL_INLINE void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
  template<time::mode TimeMode = time::mode::slew>
  VESSL_INLINE void generate(array<T> output) { proc_gen<TimeMode, false>(output, output); }

  template<time::mode TimeMode = time::mode::slew>
  VESSL_INLINE void process(array<T> input, array<T> output) { proc_gen<TimeMode, true>(input, output); }
    
protected:
  VESSL_INLINE parameter element_at(size_t index) const override
  {
    parameter p[num] = { enabled(), position(), duration(), rate() };
    return p[index];
  }

private:
  // shared routine used by templated processing and generation methods
  template<time::mode TimeMode, bool UseInput>
  void proc_gen(array<T> input, array<T> output);
    
  struct
  {
    binary_p enabled;
    analog_p position;
    duration_p duration;
    analog_p rate;
  } params_;
    
  sample::delay_line<T> delay_line_;
  analog_t phase_;
  // used to crossfade between the incoming signal and the freeze signal when enabled changes.
  math::easing::smoother<analog_t> crossfade_;
  analog_t freeze_delay_, freeze_size_;
  analog_t read_rate_, dt_;
};

// unit for use with filter types defined in the filtering namespace.
// For example, for a 2 stage biquad low pass filter use:
// filter<float, filtering::biquad<2>::lowPass>
template<typename T, template<typename> typename H>
class filter : public unit_processor<T>, protected plist<3>
{
public:
  typedef H<T> function;
    
  explicit filter(analog_t sample_rate) : unit_processor<T>(), sample_rate_(sample_rate)
  {
    params_.fhz.value = sample_rate;
    params_.q.value = 1;
    params_.emphasis.value = gain_t::from_decibels(0);
  }
    
  filter(analog_t sample_rate, analog_t freq_in_hz, analog_t kyu = filtering::q::butterworth<analog_t>()
    , gain_t emphasis = gain_t::from_decibels(0) ) 
  : unit_processor<T>(), sample_rate_(sample_rate)
  {
    params_.fhz.value = freq_in_hz;
    params_.q.value = kyu;
    params_.emphasis.value = emphasis;
  }
    
  void set_sample_rate(analog_t sample_rate) override { sample_rate_ = sample_rate; }
  const parameter_list& parameters() const override { return *this; }

  parameter fhz() const { return params_.fhz("fHz", 'f'); }
  parameter q() const { return params_.q("q", 'q'); }
  // unused by some filter types (see filtering section)
  parameter emphasis() const { return params_.emphasis("emphasis", 'e'); }
    
  T process(const T& in) override;
  void process(array<T> in, array<T> out) override;

protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = {fhz(), q(), emphasis()}; return p[index];
  }
    
private:
  struct
  {
    analog_p fhz;
    analog_p q;
    gain_p   emphasis;
  } params_;
  function func_;
  analog_t sample_rate_;
};

// designed to work with floating point types.
template<typename T, uint32_t MaxBits>
class bitcrush : public unit_processor<T>, protected plist<3>
{
public:
  bitcrush(analog_t sample_rate, analog_t bit_rate, analog_t bit_depth = MaxBits)
    : unit_processor<T>(), prev_input_(0), curr_sample_(0), rate_alpha_(0), dt_(1.0f/sample_rate)
  {
    params_.bitRate.value = bit_rate;
    params_.bitDepth.value = bit_depth;
  }
    
  void set_sample_rate(analog_t sampleRate) override { dt_ = 1.0f / sampleRate;}
  const parameter_list& parameters() const override { return *this; }

  parameter rate() const { return params_.bitRate("bit rate", 'r'); }
  parameter depth() const { return params_.bitDepth("bit depth", 'd'); }
  parameter mangle() const { return params_.mangle("mangle", 'm'); }

  T process(const T& in) override;

  using unit_processor<T>::process;
    
protected:
  parameter element_at(size_t index) const override
  {
    parameter p[num] = { rate(), depth(), mangle() };
    return p[index];
  }
    
private:
  struct
  {
    analog_p bitRate;
    analog_p bitDepth;
    binary_p mangle;
  } params_;
  T prev_input_;
  T curr_sample_;
  analog_t rate_alpha_;
  analog_t dt_;
};

} // namespace processors
} // namespace vessl

#include "vessl_generators.inl"
#include "vessl_processors.inl"
