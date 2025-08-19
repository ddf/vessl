////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2019 Damien Quartz
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
///////////////////////////////////////////////////////////////////////////////////

// ReSharper disable CppClangTidyPortabilityTemplateVirtualMemberFunction
#pragma once
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

// When built for ARM Cortex-M processor series,
// we provide template specializations that use the optimized CMSIS library:
// http://www.keil.com/pack/doc/CMSIS/General/html/index.html
#ifdef ARM_CORTEX
#include "arm_math.h" 
#endif //ARM_CORTEX

#ifndef PI
#define PI 3.14159265358979323846
#endif

// mainly to get Rider to shut up about not being able to find assert even though we include <cassert>
#ifndef NDEBUG
#ifndef assert
static void assert(bool condition) { }
#endif
#endif

#define VASSERT(cond, msg) assert((void(msg), cond))

// Note: In all classes using typename T, it is assumed to be POD and to have support for all arithmetical operators
namespace vessl
{
  template<typename T>
  class source
  {
  public:
    source() = default;
    virtual ~source() = default;
    source(const source&) = delete;
    source(const source&&) = delete;
    source& operator=(const source&) = delete;
    source& operator=(source&&) = delete;

    virtual bool isEmpty() const = 0;
    virtual T read() = 0;

    explicit operator bool() const { return !isEmpty(); }
  };

  template<typename T>
  class sink
  {
  public:
    sink() = default;
    virtual ~sink() = default;
    sink(const sink&) = delete;
    sink(const sink&&) = delete;
    sink& operator=(const sink&) = delete;
    sink& operator=(sink&&) = delete;

    virtual bool isFull() const = 0;
    virtual void write(const T& value) = 0;

    explicit operator bool() const { return !isFull(); }
  };

  template<typename T>
  class generator : public source<T>
  {
  public:
    generator() = default;
    virtual T generate() = 0;

    // by default, we assume an endless source of data
    bool isEmpty() const override { return false; }
    T read() override { return generate(); }
  };

  template<typename T>
  class array
  {
  protected:
    T* data;
    size_t size;

  public:
    array() : data(nullptr), size(0) {}
    array(T* data, size_t size) : data(data), size(size) {}

    T* getData() { return data; }
    const T* getData() const { return data; }
    size_t getSize() const { return size; }
    bool isEmpty() const { return size == 0; }
    T& operator[](size_t index) { return data[index]; }
    const T& operator[](size_t index) const { return data[index]; }

    // ranged-based for support
    T* begin() { return data; }
    T* end() { return data + size; }

    class reader final : public source<T>
    {
    protected:
      const T* begin;
      const T* head;
      const T* end;

    public:
      explicit reader(array source) : source<T>(), begin(source.data), head(source.data), end(source.data + source.size) {}
      reader(T* data, size_t size) : source<T>(), begin(data), head(data), end(data + size) {}

      // source methods
      bool isEmpty() const override { return head == end; }
      T read() override { return *head++; }

      size_t available() const { return end - head; }

      T peek() const { return *head; }
      const T* operator*() const { return head; }

      reader reset() { head = begin; return *this; }
    };

    class writer final : public sink<T>
    {
    protected:
      T* head;
      const T* end;

    public:
      explicit writer(array source) : sink<T>(), head(source.data), end(source.data + source.size) {}
      writer(T* data, size_t size) : sink<T>(), head(data), end(data + size) {}

      bool isFull() const override { return head == end; }
      void write(const T& v) override { *head++ = v; }

      size_t available() const { return end - head; }

      // block copy the entire contents of reader into this writer, returns this.
      // writer must have enough space for the contents of reader.
      writer operator<<(const reader& r);
    };

    reader getReader() const { return reader(*this); }
    writer getWriter() { return writer(*this); }

    array operator<<(array copyFrom);

    // returns dest
    array add(array other, array dest);
    // returns this
    array add(array other) { return add(other, *this); }
    // returns dest
    array scale(T value, array dest);
    // returns this
    array scale(T value) { return scale(value, *this); }
  };

  template<typename T>
  T* begin(array<T>& arr) { return arr.begin(); }

  template<typename T>
  T* end(array<T>& arr) { return arr.end(); }

  template<typename T>
  class processor
  {
  public:
    processor() = default;
    virtual ~processor() = default;
    processor(const processor&) = delete;
    processor(const processor&&) = delete;
    processor& operator=(const processor&) = delete;
    processor& operator=(processor&&) = delete;

    virtual T process(const T& in) = 0;
    virtual void process(source<T>& in, sink<T>& out);
    virtual void process(array<T> in, array<T> out);
  };

  template<typename T>
  class ring : array<T>
  {
    T* head;

  public:
    ring(T* inData, size_t inSize) : array<T>(inData, inSize), head(inData + inSize - 1) {}

    // expose direct access to underlying array data
    using array<T>::getData;
    using array<T>::getSize;

    void write(const T& v);
    size_t getWriteIndex() const { return head - array<T>::data; }

    ring operator<<(typename array<T>::reader r);
  };

  class parameter
  {
  public:
    enum class type : uint8_t
    {
      binary = 0, // on/off
      digital = 1, // integral values
      analog = 2, // floating point values
      // space for more built-ins

      user = UINT8_MAX // user provided type, stored as a void*
    };

    parameter(const char* name, type type);
    parameter(const char* name, bool value) : pn(name), pt(type::binary) { pv.b = value; }
    parameter(const char* name, long value) : pn(name), pt(type::digital) { pv.i = value; }
    parameter(const char* name, double value) : pn(name), pt(type::analog) { pv.a = value; }
    parameter(const char* name, void* value) : pn(name), pt(type::user) { pv.u = value; }

    const char* getName() const { return pn; }
    type getType() const { return pt; }

    // for the sane, static cast from our stored type to T
    template<typename T>
    T read() const;

    // read the state of a binary parameter, including sample delay
    bool read(uint32_t* sampleDelay) const;

    template<typename T>
    parameter& write(const T& value);

    // set the state of a binary parameter with a sample delay
    parameter& write(bool gate, uint32_t sampleDelay);

    // overloading dereference with float conversion because it will be used so often
    float operator*() const { return read<float>(); }

    template<typename T>
    T operator!() const;

    template<typename T>
    T operator-() const;

  private:
    const char* pn;

    union
    {
      size_t  b;
      int64_t i;
      double  a;
      void*   u;
    } pv;

    type pt;
  };

  template<typename T>
  parameter& operator<<(parameter& p, const T& v) { return p.write(v); }
  
  template<typename T>
  T& operator>>(parameter& p, T& v) { v = p.read<T>(); return v; }
  
  template<typename T>
  bool operator>(const parameter& p, const T& v) { return p.read<T>() > v; }
  
  template<typename T>
  bool operator<(const parameter& p, const T& v) { return p.read<T>() < v; }
  
  template<typename T>
  bool operator==(const parameter& p, const T& v) { return p.read<T>() == v; }
  
  template<typename T>
  bool operator!=(const parameter& p, const T& v) { return p.read<T>() != v; }
  
  template<typename T>
  T operator+(const parameter& p, const T& v) { return p.read<T>() + v; }
  
  template<typename T>
  T operator+(const T& v, const parameter& p) { return p.read<T>() + v; }
  
  template<typename T>
  T operator*(const parameter& p, const T& v) { return p.read<T>() * v; }
  
  template<typename T>
  T operator*(const T& v, const parameter& p) { return p.read<T>() * v; }
  
  template<typename T>
  T operator/(const parameter& p, const T& v) { return p.read<T>() / v; }
  
  template<typename T>
  T operator/(const T& v, const parameter& p) { return v / p.read<T>(); }

  class unit : array<parameter>
  {
    const char* unitName;
    float sampleRate;
    float deltaTime;
    
  protected:
    template<size_t N>
    struct init
    {
      const char* name;
      parameter params[N];
    };

    template<size_t N>
    explicit unit(init<N>& init, float sampleRate = 1) : array(init.params, N), unitName(init.name) { setSampleRate(sampleRate); }
    unit(const char* name, parameter* params, int paramsCount, float sampleRate = 1) : array(params, paramsCount), unitName(name) { setSampleRate(sampleRate); }

  public:
    virtual ~unit() = default;
    unit(const unit&) = default;
    unit(unit&&) = default;
    unit& operator=(const unit&) = default;
    unit& operator=(unit&&) = default;

    const char* name() const { return unitName; }

    float getSampleRate() const { return sampleRate; }

    void setSampleRate(float sr);

    parameter& getParameter(size_t index) { return this->operator[](index); }
    const parameter& getParameter(size_t index) const { return this->operator[](index); }
    size_t getParameterCount() const { return getSize(); }

    using array::begin;
    using array::end;

  protected:
    virtual void onSampleRateChanged() {}
    float dt() const { return deltaTime; }
  };

  template<typename T>
  class unitGenerator : public unit, public generator<T>
  {
  protected:
    template<size_t N>
    explicit unitGenerator(init<N>& init, float sampleRate = 1) : unit(init, sampleRate), generator<T>() {}
    unitGenerator(const char* name, parameter* params, int paramsCount, float sampleRate = 1) : unit(name, params, paramsCount, sampleRate), generator<T>() {}
  };

  template<typename T>
  class unitProcessor : public unit, public processor<T>
  {
  protected:
    template<size_t N>
    explicit unitProcessor(init<N>& init, float sampleRate = 1) : unit(init, sampleRate), processor<T>() {}
    unitProcessor(const char* name, parameter* params, int paramsCount, float sampleRate = 1) : unit(name, params, paramsCount, sampleRate), processor<T>() {}
  };

  namespace interpolation
  {
    template<typename T>
    struct nearest
    {
      T operator()(const T* buffer, float fracIdx);
    };

    template<typename T>
    struct linear
    {
      T operator()(const T* buffer, float fracIdx);
    };

    template<typename T>
    struct cubic
    {
      T operator()(const T* buffer, float fracIdx);
    };
    
    template<typename T, typename I = linear<T>>
    T sample(const T* buffer, float fracIdx)
    {
      VASSERT(fracIdx >= 0, "fracIdx argument to sample must be non-negative");
      static I interpolator;
      return interpolator(buffer, fracIdx);
    }
  };

  // @todo add all functions used by vessl and call into here.
  // for the most part these should default to standard math calls,
  // so that if those calls are redefined by a platform, we can take advantage of that.
  // e.g. by including basicmaths.h before including vessl on OWL.
  namespace math
  {
    template<typename T>
    struct epsilon { static constexpr T value = std::numeric_limits<T>::epsilon(); };
    
    template<typename T>
    T max(T a, T b) { return a > b ? a : b; }

    template<typename T>
    T mod(T v, T* i) { return modf(v, i); }
    
    template<typename T>
    T wrap(T val, T low, T high)
    {
      T diff = high - low;
      while (val < low)
      {
        val += diff;
      }
      while (val > high)
      {
        val -= diff;
      }
      return val;
    }

    template<typename T>
    T wrap01(T val) { return wrap(val, T(0), T(1)); }
    
    template<typename T>
    T constrain(T val, T low, T high)
    {
      return val < low ? low : val > high ? high : val;
    }

    template<typename T>
    bool isNan(T n) { return isnan(n); }
  }

  namespace random
  {
    static constexpr int I32_MAX = RAND_MAX;
    inline void si32(unsigned seed) { srand(seed); }
    inline int i32() { return rand(); }
  }
  
  namespace easing
  {
    struct linear
    {
      float operator()(float t) const { return t; }
    };

    template<typename T, typename E = linear>
    T interp(T begin, T end, float t)
    {
      static E ease;
      return (end-begin) * ease(t) + begin;
    }
  }
  
  struct time
  {
    // convert bpm to frequency in Hz
    static constexpr float B2F = 1.0f / 60.0f;
    static constexpr float F2B = 60.f;
    
    float sampleRate;
    size_t period; // in samples

    static time fromBpm(float bpm, float sr) { return { sr, static_cast<size_t>(sr/(bpm*B2F)) }; }
    float toBpm() const { return F2B*(sampleRate/static_cast<float>(period)); }
  };
  
  // support for tempo detection of a clock signal (i.e. pulse train)
  class clockable
  {
  protected:
    time   tempo;
    size_t periodMin;
    size_t periodMax;
    size_t ticks;

    clockable(float sampleRate, size_t samplePeriodMin, size_t samplePeriodMax, float bpm = 60)
    : tempo(time::fromBpm(bpm, sampleRate)), periodMin(samplePeriodMin), periodMax(samplePeriodMax), ticks(0)
    {}

    // subclasses should call tick for every sample generated/processed
    void tick() { ++ticks; }
    void tick(size_t t) { ticks += t; }

    // subclass can override this to be notified every time they receive a clock pulse
    virtual void tock(size_t sampleDelay) {}
    
  public:
    virtual ~clockable() = default;
    clockable(const clockable&) = default;
    clockable(clockable&&) = default;
    clockable& operator=(const clockable&) = default;
    clockable& operator=(clockable&&) = default
    ;
    // users should call clock at the beginning of every clock pulse
    void clock() { tempo.period = math::constrain(ticks, periodMin, periodMax); ticks = 0; tock(0); }
    void clock(size_t sampleDelay) { tempo.period = math::constrain(ticks + sampleDelay, periodMin, periodMax); ticks = 0; tock(sampleDelay); }

    float getBpm() const { return tempo.toBpm(); }
  };
  
  // a waveform that can be evaluated using a normalized phase value
  // implementors should accept negative phase, as well as phase values outside [-1,1]
  template<typename T>
  class waveform
  {
  public:
    waveform() = default;
    virtual ~waveform() = default;
    waveform(const waveform&) = delete;
    waveform(const waveform&&) = delete;
    waveform& operator=(const waveform&) = delete;
    waveform& operator=(waveform&&) = delete;

    virtual T evaluate(float phase) const = 0;  // NOLINT(portability-template-virtual-member-function)
  };

  namespace waves
  {
    template<typename T = double>
    class sine final : public waveform<T>
    {
    public:
      sine() = default;
      T evaluate(float phase) const override { return sin(T(2 * PI * phase)); }
    };
  }

  // a fixed-sized buffer that supports sampling it a normalized phase.
  // both positive and negative phases are supported
  template<typename T, size_t N, typename I = interpolation::linear<T>>
  class wavetable final : public waveform<T>
  {
    T buffer[N + 3] = {};

  public:
    wavetable(): waveform<T>() {}

    // assumes source can provide at least N samples
    explicit wavetable(source<T>& source);
    explicit wavetable(const waveform<T>& waveform);

    // ReSharper disable once CppMemberFunctionMayBeStatic
    size_t size() const { return N; }

    T get(const size_t i) const { return buffer[i + 1]; }
    void set(const size_t i, T val);
    
    // implement waveform:
    T evaluate(float phase) const override;
  };

  // @todo refactor this to work like interpolation?
  enum class noiseTint : uint8_t
  {
    white,
    pink,
    red,
    brown,
  };

  template<typename T>
  class noiseGenerator final : public unitGenerator<T>
  {
    unit::init<2> init = {
      "noise",
      {
        parameter("tint", parameter::type::digital),
        parameter("rate", parameter::type::analog)
      }
    };
    float step;
    T nz[2];

    using unit::dt;

  public:
    explicit noiseGenerator(float sampleRate, noiseTint withTint = noiseTint::white)
    : unitGenerator<T>(init, sampleRate), step(0)
    {
      nz[0] = nz[1] = next(withTint, 0);
      tint() << withTint;
      rate() << 1.0;
    }
    noiseGenerator(const noiseGenerator&) = default;
    noiseGenerator(noiseGenerator&&) = default;
    noiseGenerator& operator=(const noiseGenerator&) = default;
    noiseGenerator& operator=(noiseGenerator&&) = default;
    ~noiseGenerator() override = default;

    parameter& tint() { return init.params[0]; }
    parameter& rate() { return init.params[1]; }

    static T next(noiseTint tint, double dt);

    T generate() override;
  };

  // unit that generates a linear ramp from one value to another over a duration of seconds
  // @todo implement easings above and add that as a template parameter
  template<typename T>
  class ramp final : public unitGenerator<T>
  {
    T mFrom, mTo;
    unit::init<4> init = {
      "ramp",
      {
        parameter("from", &mFrom),
        parameter("to", &mTo),
        parameter("duration", parameter::type::analog),
        parameter("eor", parameter::type::binary)
      }
    };
    float t;

    using unit::dt;

  public:
    explicit ramp(float sampleRate, float durationInSeconds = 0, T fromValue = T(0), T toValue = T(0))
    : unitGenerator<T>(init, sampleRate), t(0)
    {
      // this is fun
      from() << fromValue;
      to() << toValue;
      duration() << durationInSeconds;
    }
    ramp(const ramp&) = default;
    ramp(ramp&&) = default;
    ramp& operator=(const ramp&) = default;
    ramp& operator=(ramp&&) = default;
    ~ramp() override = default;

    // ins
    parameter& from() { return init.params[0]; }
    parameter& to() { return init.params[1]; }
    parameter& duration() { return init.params[2]; }

    // outs
    const parameter& eor() const { return init.params[3]; }
    // could also add t as an out.

    bool isActive() const { return !eor().template read<bool>(); }
    void trigger();
    T generate() override;

  private:
    parameter& eorw() { return init.params[3]; }
  };

  template<typename T>
  class oscil final : public unitGenerator<T>
  {
    unit::init<4> init = {
      "oscil", {
        parameter("frequency", parameter::type::analog),
        parameter("fm (lin)", parameter::type::analog),
        parameter("fm (v/oct)", parameter::type::analog),
        parameter("phase mod", parameter::type::analog)
      }
    };

    waveform<T>* wave;
    float phase;

    using unit::dt;

  public:
    explicit oscil(float sampleRate, waveform<T>& wave, double freqInHz = 440)
    : unitGenerator<T>(init, sampleRate), wave(&wave), phase(0)
    {
      fHz() << freqInHz;
    }
    ~oscil() override = default;
    oscil(const oscil&) = default;
    oscil(oscil&&) = default;
    oscil& operator=(const oscil&) = default;
    oscil& operator=(oscil&&) = default;

    // frequency in Hz without FM applied
    parameter& fHz() { return init.params[0]; }
    // linear frequency modulation
    parameter& fmLin() { return init.params[1]; }
    // v/oct (exponential) frequency modulation
    parameter& fmExp() { return init.params[2]; }
    // phase modulation
    parameter& pm() { return init.params[3]; }

    T generate() override;

    void reset() { phase = 0; }
  };

  template<typename T>
  class delayline final : public ring<T>, public waveform<T>
  {
  public:
    delayline(T* inData, size_t inSize) : ring<T>(inData, inSize) {}

    using ring<T>::getData;
    using ring<T>::getSize;
    using ring<T>::getWriteIndex;

    // reads behind the write head with sampleDelay (i.e. the ith sample previously written)
    // where a delay of 0 samples will give the most recently written value.
    T read(size_t sampleDelay);

    // reads behind the write head with a fractional sampleDelay and given interpolation
    template<typename I = interpolation::linear<T>>
    T read(float sampleDelay) const;

    // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
    T evaluate(float phase) const override;
  };

  template<typename T, typename I = interpolation::linear<T>>
  class delay : public unitProcessor<T>
  {
    unit::init<2> init = {
      "delay",
      {
        parameter("time", parameter::type::analog),
        parameter("feedback", parameter::type::analog)
      }
    };
    delayline<T> buffer;
    float mTime;

    using unit::dt;

  public:
    delay(array<T> buffer, float sampleRate, float delayInSeconds = 0, float feedbackAmount = 0)
    : unitProcessor<T>(init, sampleRate), buffer(buffer.getData(), buffer.getSize()), mTime(delayInSeconds)
    {
      time() << delayInSeconds;
      feedback() << feedbackAmount;
    }

    /// amount of delay time in seconds
    parameter& time() { return init.params[0]; }
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    parameter& feedback() { return init.params[1]; }

    T process(const T& in) override;

    using processor<T>::process;
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // implementation
  //

  template<typename T>
  void processor<T>::process(source<T>& in, sink<T>& out)
  {
    while (in && out)
    {
      out.write(process(in.read()));
    }
  }

  template<typename T>
  void processor<T>::process(array<T> in, array<T> out)
  {
    auto r = in.getReader();
    auto w = out.getWriter();
    process(r, w);
  }

  template<typename T>
  typename array<T>::writer array<T>::writer::operator<<(const reader& r)
  {
    size_t rsz = r.available();
    size_t wsz = available();
    VASSERT(wsz >= rsz, "Not enough space in writer for the contents of reader");
    const T* rh = *r;
    memcpy(static_cast<void*>(head), static_cast<const void*>(rh), rsz * sizeof(T));
    head += rsz;
    return writer(head, wsz - rsz);
  }
  
  template<typename T>
  array<T> array<T>::operator<<(array copyFrom)
  {
    writer w(*this);
    w << reader(copyFrom);
    return *this;
  }

  template<typename T>
  array<T> array<T>::add(array other, array dest)
  {
    VASSERT(size == other.size && size <= dest.size, "arrays are have different lengths or destination is too small");
    reader a(data, size);
    reader b(other.data, other.size);
    writer c(dest.data, dest.size);
    while (a)
    {
      c << a.read() + b.read();
    }
    return dest;
  }

  template<typename T>
  array<T> array<T>::scale(T value, array dest)
  {
    VASSERT(size <= dest.size, "destination size is too small");
    reader a(data, size);
    writer b(dest.data, dest.size);
    while (a)
    {
      b << a.read() * value;
    }
    return dest;
  }

  template<typename T>
  typename array<T>::writer& operator<<(typename array<T>::writer& w, const T& v)
  {
    w.write(v);
    return w;
  }
  
  template<typename T>
  void ring<T>::write(const T& v)
  {
    *head++ = v;
    if (head == array<T>::end())
    {
      head = array<T>::begin();
    }
  }

  template<typename T>
  ring<T> ring<T>::operator<<(typename array<T>::reader r)
  {
    VASSERT(r.available() < getSize(), "reader size is larger than ring size");
    while (r)
    {
      write(r.read());
    }
    return *this;
  }
  
  inline parameter::parameter(const char* name, type type): pn(name), pt(type)
  {
    switch (pt)
    {
      case type::binary: pv.b = 0; break;
      case type::digital: pv.i = 0; break;
      case type::analog: pv.a = 0; break;
      case type::user: pv.u = nullptr; break;
    }
  }

  inline bool parameter::read(uint32_t* sampleDelay) const
  {
    *sampleDelay = static_cast<uint32_t>(pv.b >> 1);
    return pv.b & 1;
  }

  inline parameter& parameter::write(bool gate, uint32_t sampleDelay)
  {
    pv.b = gate | sampleDelay << 1;
    return *this;
  }

  template<typename T>
  T parameter::read() const
  {
    // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
    switch (pt)
    {
      case type::binary:
      {
        // if we were set with a sample delay, decrement the sample delay
        // and report the opposite value of our first bit
        size_t state = pv.b & 1;
        uint32_t sampleDelay = static_cast<uint32_t>(pv.b >> 1);
        if (sampleDelay == 0) { return static_cast<T>(state); }
        const_cast<parameter*>(this)->pv.b = state | ((sampleDelay - 1) << 1);
        return static_cast<T>(!state);
      }
      case type::digital: return static_cast<T>(pv.i);
      case type::analog: return static_cast<T>(pv.a);
      case type::user: return pv.u ? *static_cast<T*>(pv.u) : T();
    }
    return T();
  }

  template<typename T>
  parameter& parameter::write(const T& value)
  {
    // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
    switch (pt)
    {
      case type::binary: pv.b = static_cast<bool>(value); break;
      case type::digital: pv.i = static_cast<int64_t>(value); break;
      case type::analog: pv.a = static_cast<double>(value); break;
      case type::user: if (pv.u) { *static_cast<T*>(pv.u) = value; } break;
    }
    return *this;
  }

  template<typename T>
  T parameter::operator!() const
  {
    // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
    switch (pt)
    {
      case type::binary: return !read<bool>();
      case type::digital: return !pv.i;
      case type::analog: return !pv.a;
      case type::user: return !(pv.u ? *pv.u : T());
    }
    return !T();
  }

  template<typename T>
  T parameter::operator-() const
  {
    // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
    switch (pt)
    {
      case type::binary: return !read<bool>();
      case type::digital: return -pv.i;
      case type::analog: return -pv.a;
      case type::user: return -(pv.u ? *pv.u : T());
    }
    return -T();
  }
  
  inline void unit::setSampleRate(float sr)
  {
    sampleRate = sr;
    deltaTime = 1.0f / sr;
    onSampleRateChanged();
  }

  template<typename T>
  T interpolation::nearest<T>::operator()(const T* buffer, float fracIdx)
  {
    return buffer[static_cast<size_t>(round(fracIdx))];
  }

  template<typename T>
  T interpolation::linear<T>::operator()(const T* buffer, float fracIdx)
  {
    float idx;
    float frac = math::mod(fracIdx, &idx);
    size_t x0 = static_cast<size_t>(idx);
    return buffer[x0] + (buffer[x0 + 1] - buffer[x0]) * frac;
  }

  template<typename T>
  T interpolation::cubic<T>::operator()(const T* buffer, float fracIdx)
  {
    static const T DIV6 = static_cast<T>(1. / 6.);
    static const T DIV2 = static_cast<T>(0.5);

    float idx;
    float f = math::mod(fracIdx, &idx);
    float fm1 = f - 1.f;
    float fm2 = f - 2.f;
    float fp1 = f + 1.f;
    size_t x0 = idx;
    return -f * fm1 * fm2 * DIV6 * buffer[x0 - 1] + fp1 * fm1 * fm2 * DIV2 * buffer[x0] - fp1 * f * fm2 * DIV2 * buffer[x0 + 1] + fp1 * f * fm1 * DIV6 * buffer[x0 + 2];
  }

  template<typename T, size_t N, typename I>
  wavetable<T, N, I>::wavetable(source<T>& source): waveform<T>()
  {
    for (int i = 0; i < N; ++i)
    {
      buffer[i + 1] = source.read();
    }

    // configure extra values on the ends of the buffer
    // so that sampling buffers that are periodic waveforms will work correctly
    buffer[0] = buffer[N];
    buffer[N + 1] = buffer[1];
    buffer[N + 2] = buffer[2];
  }

  template<typename T, size_t N, typename I>
  wavetable<T, N, I>::wavetable(const waveform<T>& waveform): waveform<T>()
  {
    float phase = 0;
    float step = 1.0f / N;
    for (size_t i = 0; i < N; ++i)
    {
      buffer[i + 1] = waveform.evaluate(phase);
      phase += step;
    }

    // configure extra values on the ends of the buffer
    // so that sampling buffers that are periodic waveforms will work correctly
    buffer[0] = buffer[N];
    buffer[N + 1] = buffer[1];
    buffer[N + 2] = buffer[2];
  }

  template<typename T, size_t N, typename I>
  void wavetable<T, N, I>::set(const size_t i, T val)
  {
    buffer[i + 1] = val;
    // also update wrap-around values
    switch (i)
    {
      case 0: buffer[N + 1] = buffer[1];
        break;
      case 1: buffer[N + 2] = buffer[2];
        break;
      case N - 1: buffer[0] = buffer[N];
        break;
      default: break;
    }
  }

  template<typename T, size_t N, typename I>
  T wavetable<T, N, I>::evaluate(float phase) const
  {
    float idx = phase * N;
    while (idx > N) { idx -= N; }
    while (idx < 0) { idx += N; }
    return interpolation::sample<T, I>(buffer, idx);
  }

  template<typename T>
  T noiseGenerator<T>::next(noiseTint tint, double dt)
  {
    switch (tint)
    {
      case noiseTint::white:
      {
        return 2 * (static_cast<T>(random::i32()) / random::I32_MAX) - 1;
      }

      // #TODO: this contains clicks when run at audio frequency and I'm not sure why.
      // Need to either read up on how to do this properly and write a new implementation,
      // or find some open source code that can be used without causing license issues.
      case noiseTint::red:
      case noiseTint::brown:
      {
        static const T RC = static_cast<T>(1) / (PI * 200);
        static const T AC = 6.2;
        static T prev = 0;
        T alpha = dt / (dt + RC);
        T white = 2 * (static_cast<T>(random::i32()) / random::I32_MAX) - 1;
        prev = easing::interp(prev, white, alpha);
        return prev * AC;
      }

      // This is the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
      // Would be good to dig into the improvements on the algorithm mentioned later in the article.
      case noiseTint::pink:
      {
        static constexpr int RANGE = 128;
        static constexpr int MAX_KEY = 0x1f;
        static int key = 0;
        static T maxSum = 90;
        static int whiteValues[6] = {
          random::i32() % (RANGE / 6), random::i32() % (RANGE / 6), random::i32() % (RANGE / 6),
          random::i32() % (RANGE / 6), random::i32() % (RANGE / 6), random::i32() % (RANGE / 6) };

        int lastKey = key;
        T sum = 0;
        key = key == MAX_KEY ? 0 : ++key;

        int diff = lastKey ^ key;
        for (int i = 0; i < 6; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            whiteValues[i] = random::i32() % (RANGE / 6);
          }
          sum += whiteValues[i];
        }
        maxSum = math::max(sum, maxSum);
        T n = 2 * (sum / maxSum) - 1;
        vassert(!math::isNan(n) && "pink noise generated nan");
        return n;
      }
    }

    return 0;
  }

  template<typename T>
  T noiseGenerator<T>::generate()
  {
    step += dt() * rate().read();
    float alpha = step / dt();
    if (alpha >= 1)
    {
      noiseTint tnt = static_cast<noiseTint>(tint());
      nz[0] = nz[1];
      nz[1] = next(tnt, dt());
      alpha = math::wrap01(alpha);
      step = alpha * dt();
    }
    return lerp(nz[0], nz[1], alpha);
  }

  template<typename T>
  void ramp<T>::trigger()
  {
    if (duration() > math::epsilon<float>::value)
    {
      t = 0;
      eorw() << T(0);
    }
    else
    {
      t = 1;
      eorw() << T(1);
    }
  }

  template<typename T>
  T ramp<T>::generate()
  {
    float lt = t;
    if (isActive())
    {
      float dinv = 1.0f / duration();
      t += dt() * dinv;
      if (t >= 1)
      {
        eorw() << T(1);
      }
    }
    return easing::interp(mFrom, mTo, lt);
  }

  template<typename T>
  T oscil<T>::generate()
  {
    T val = wave->evaluate(phase + pm());
    phase += (fHz() * exp2(*fmExp()) + fmLin()) * dt();
    phase = math::wrap01(phase);
    return val;
  }

  template<typename T>
  T delayline<T>::read(size_t sampleDelay)
  {
    assert(sampleDelay < getSize() - 1);
    size_t idx = getWriteIndex() + 1 + sampleDelay;
    return getData[idx % getSize()];
  }

  template<typename T>
  template<typename I>
  T delayline<T>::read(float sampleDelay) const
  {
    float fSize = static_cast<float>(getSize());
    sampleDelay = fSize - sampleDelay;
    assert(sampleDelay >= 0 && sampleDelay < fSize - 1);
    sampleDelay = static_cast<float>(getWriteIndex() + 1) + sampleDelay;
    float idx;
    float f = math::mod(sampleDelay, &idx);
    size_t x0 = static_cast<size_t>(idx) % getSize();
    size_t x1 = (x0 + 1) % getSize();
    size_t x2 = (x0 + 2) % getSize();
    const T* data = getData();
    T s[3] = { data[x0], data[x1], data[x2] };
    return interpolation::sample<T, I>(s, f);
  }

  template<typename T>
  T delayline<T>::evaluate(float phase) const
  {
    float fSize = static_cast<float>(getSize());
    phase = math::wrap(phase, -1.f, 1.f);
    float sampleDelay = phase > 0 ? (1.0f - phase) * fSize : -phase * fSize;
    return read(sampleDelay);
  }

  template<typename T, typename I>
  T delay<T, I>::process(const T& in)
  {
    // smooth time parameter to prevent crunchiness when it is noisy or changes by large amounts
    mTime = easing::interp(mTime, *time(), dt() * 10);
    // delay time in samples
    float dts = math::constrain(mTime * unit::getSampleRate(), 0.f, static_cast<float>(buffer.getSize()-1));
    float fbk = math::constrain(feedback() >> fbk, -1.f, 1.f);
    T s = buffer.template read<I>(dts);
    buffer.write(in + s * fbk);
    return s;
  }

  template<>
  inline float math::wrap01(float v)
  {
    float i;
    return modf(v, &i);
  }

#ifdef ARM_CORTEX
  template<>
  array<float> array<float>::add(array other, array dest)
  {
    assert(size == other.size && size <= dest.size);
    arm_add_f32(data, other.data, dest.data, size);
    return dest;
  }

  template<>
  array<float> array<float>::scale(float value, array dest)
  {
    assert(size <= dest.size);
    arm_scale_f32(data, value, dest.data, size);
    return dest;
  }
#endif
}
