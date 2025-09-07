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
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

// When built for ARM Cortex-M processor series,
// we provide template specializations that use the optimized CMSIS library:
// http://www.keil.com/pack/doc/CMSIS/General/html/index.html
#ifdef ARM_CORTEX
#include "arm_math.h" 
#endif //ARM_CORTEX

// because some people like to redefine these math functions with macros
#ifdef sqrt
#undef sqrt
#endif

#ifdef pow
#undef pow
#endif

#ifdef round
#undef round
#endif

#ifdef sin
#undef sin
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
  // @todo use these throughout
  using binary_t = bool;
  using integral_t = int64_t;
  using analog_t = double;
  
  template<typename T>
  class source
  {
  public:
    source() = default;
    virtual ~source() = default;
    source(const source&) = default;
    source(source&&) = default;
    source& operator=(const source&) = default;
    source& operator=(source&&) = default;

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
    sink(const sink&) = default;
    sink(sink&&) = default;
    sink& operator=(const sink&) = default;
    sink& operator=(sink&&) = default;

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
      reader(const T* data, size_t size) : source<T>(), begin(data), head(data), end(data + size) {}

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

    void fill(T value);
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
  
  namespace frame
  {
    template<typename T, size_t N>
    struct channels : array<T>
    {
      T samples[N];

      channels() : array<T>(samples, N) { array<T>::fill(0); }
    };

    template<typename T>
    class mono : public channels<T, 1>
    {
      using channels<T,1>::samples;
  
    public:
      T& value() { return samples[0]; }
      const T& value() const { return samples[0]; }
    };

    template<typename T>
    class stereo : public channels<T, 2>
    {
      using channels<T, 2>::samples;
      
    public:
      stereo() {}
      stereo(T l, T r) { samples[0] = l, samples[1] = r; }
    
      T& left() { return samples[0]; }
      const T& left() const { return samples[0]; }
      T& right() { return samples[1]; }
      const T& right() const { return samples[1]; }
    };
  };

  template<typename T>
  class processor
  {
  public:
    processor() = default;
    virtual ~processor() = default;
    processor(const processor&) = delete;
    processor(processor&&) = delete;
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
    void setWriteIndex(size_t index) { head = array<T>::data + index%array<T>::size; }

    ring operator<<(typename array<T>::reader r);
  };

  class parameter
  {
  public:
    enum class type : uint8_t
    {
      binary = 0, // on/off (bool)
      digital = 1, // integral values (int64_t)
      analog = 2, // floating point values (double)
      // space for more built-ins

      // user provided type, stored as a void*, must be convertible to all other parameter types.
      // even if the conversion is meaningless.
      user = UINT8_MAX 
    };

    parameter(const char* name, type type);
    parameter(const char* name, void* userData);

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

    parameter& operator<<(const parameter& rhs)
    {
      switch (rhs.pt)
      {
        case type::binary: this->write(rhs.pv.b); break;
        case type::digital: this->write(rhs.pv.i); break;
        case type::analog: this->write(rhs.pv.a); break;
        case type::user: VASSERT(false, "Can't copy a user parameter without knowing the type!"); break;
      }

      return *this;
    }

    // overloading dereference with float conversion because it will be used so often
    float operator*() const { return read<float>(); }
    explicit operator bool() const { return read<bool>(); }

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
  
  // for the most part these resolve to standard math calls,
  // template specializations are provided for ARM_CORTEX where applicable.
  namespace math
  {
    template<class T>
    constexpr T e() { return static_cast<T>(2.71828182845904523536); }
    
    template<class T>
    constexpr T pi() { return static_cast<T>(3.1415926535897932385); }

    template<class T>
    constexpr T twoPi() { return pi<T>() * 2; }

    template<typename T>
    T abs(T val)
    {
      return ::abs(val);
    }
    
    template<typename T>
    T constrain(T val, T low, T high)
    {
      return val < low ? low : val > high ? high : val;
    }
    
    template<typename T>
    T epsilon() { return std::numeric_limits<T>::epsilon(); }

    template<typename T>
    T exp10(T v) { return ::exp10(v); }

    template<typename T>
    T log10(T v) { return ::log10(v); }
    
    template<typename T>
    T max(T a, T b) { return a > b ? a : b; }

    template<typename T>
    T mod(T v, T* i) { return ::modf(v, i); }

    template<typename T>
    T pow(T x, T y) { return ::pow(x, y); }

    template<typename T>
    T round(T x) { return ::round(x); }

    template<typename T>
    T sin(T x) { return ::sin(x); }

    template<typename T>
    T sqrt(T x) { return ::sqrt(x); }

    template<typename T>
    T tan(T x) { return ::tan(x); }
    
    template<typename T>
    T wrap(T val, T low, T high)
    {
      // @todo probably a way to do this without while loops.
      T diff = high - low;
      while (val < low) { val += diff; }
      while (val > high) { val -= diff; }
      return val;
    }

    // @todo use modf here
    template<typename T>
    T wrap01(T val) { return wrap(val, T(0), T(1)); }

    template<typename T>
    bool isNan(T n) { return isnan(n); }
  }

  namespace random
  {
    static constexpr int I32_MAX = RAND_MAX;
    inline void si32(unsigned seed) { srand(seed); }
    inline int i32() { return rand(); }
  }

  template<typename T>
  class gain
  {
    T value;
  public:
    static T decibelsToScale(T db) { return math::exp10(db*0.05);}
    static T scaleToDecibels(T scale) { return math::log10(scale)*20.0;}
    static gain fromScale(float scale) { return {scale}; }
    static gain fromDecibels(float dB) { return {decibelsToScale(dB)}; }

    T toScale() const { return value; }
    T toDecibels() const { return scaleToDecibels(value); }
  };
  
  namespace easing
  {
    struct linear
    {
      float operator()(float t) const { return t; }
    };

    namespace quad
    {
      struct in { float operator()(float t) const { return t * t; } };
      struct out { float operator()(float t) const { return 1.0f - (1.0f - t) * (1.0f - t); } };
      struct inOut { float operator()(float t) const { return t < 0.5f ? 2*t*t : 1.0f - math::pow(-2.0f*t+2.0f, 2.0f) * 0.5f;  } };
      struct outIn { static out qo; float operator()(float t) const { return t < 0.5f ? qo(2.0f*t) * 0.5f : 1.0f - qo(2.0f*t) * 0.5f; } };
    }

    namespace expo
    {
      struct in { float operator()(float t) const { return t <= math::epsilon<float>() ? 0.f : math::pow(2.f, 10*t-10); } };
      struct out { float operator()(float t) const { return t >= 1.0f - math::epsilon<float>() ? 1.0f : 1.0f - math::pow(2.0f, -10*t);} };
      struct inOut
      {
        float operator()(float t) const
        {
          return t <= math::epsilon<float>() ? 0.0f
          : t >= 1.0f - math::epsilon<float>() ? 1.0f
          : t < 0.5f ? math::pow(2.f, 20*t-10)*0.5f
          : 2.0f - math::pow(2.f, -20*t+10)*0.5f;
        }
      };
    }

    // easing is first parameter so T can be deduced
    template<typename E, typename T>
    T interp(T begin, T end, float t)
    {
      static E ease;
      return (end-begin) * ease(t) + begin;
    }

    template<typename T>
    T lerp(T begin, T end, float t) { return interp<linear, T>(begin, end, t); }
  }

  namespace filtering
  {
    template<typename T>
    struct data
    {
      array<T> coeff;
      array<T> state;

      data(T* coeffData, size_t coeffSize, T* stateData, size_t stateSize)
      : coeff(coeffData, coeffSize), state(stateData, stateSize)
      { coeff.fill(0); state.fill(0); }
      
      size_t getCoeffSize() const { return coeff.getSize(); }
      size_t getStateSize() const { return state.getSize(); }
    };

    // based on https://www.earlevel.com/main/2012/11/26/biquad-c-source-code/ 
    namespace biquad
    {
      template<size_t STAGES>
      struct type
      {
        static constexpr size_t stages = STAGES;
        
        template<typename T>
        void copy(T* coeff)
        {
          size_t cosz = 5;
          array<T> src(coeff, cosz);
          for (int i = 1; i < stages; i++)
          {
            array<T> dst(coeff + cosz*i, cosz);
            dst << src;
          }
        }
      };

      template<size_t STAGES = 1>
      struct lp : type<STAGES>
      {
        template<typename T>
        void operator()(T* coeff, T omega, T q, size_t stages = 1)
        {
          T K = math::tan(omega);
          T ralpha = 1 / (1 + K / q + K * K);
          coeff[0] = K * K * ralpha;
          coeff[1] = 2 * coeff[0];
          coeff[2] = coeff[0];
          coeff[3] = - 2 * (K * K - 1) * ralpha;
          coeff[4] = - (1 - K / q + K * K) * ralpha;
          // copy calculated coefficients to all stages
          if (stages > 1)
          {
            type<STAGES>::copy(coeff);
          }
        }
      };

      template<size_t STAGES = 1>
      struct hp : type<STAGES>
      { 
        template<typename T>
        void operator()(T* coeff, T omega, T q)
        {
          T K = math::tan(omega);
          T ralpha = 1 / (1 + K / q + K * K);
          coeff[0] = 1 * ralpha;
          coeff[1] = -2 * coeff[0];
          coeff[2] = coeff[0];
          coeff[3] = - 2 * (K * K - 1) * ralpha;
          coeff[4] = - (1 - K / q + K * K) * ralpha;
          // copy calculated coefficients to all stages
          if (STAGES > 1)
          {
            type<STAGES>::copy(coeff);
          }
        }
      };

      template<typename T, size_t STATES, size_t STAGES>
      struct cascade : data<T>
      {
      private:
        T co[5*STAGES];
        T st[STATES*STAGES];
      public:
        cascade() : data<T>(co, 5*STAGES, st, STATES*STAGES) {}
      };

      // @todo want to be able to have class... Fs in the template parameter list
      // so that it is possible to declare df2T types that mix biquad::types
      // e.g.: df2T<float, bp<2>, notch<2>>
      // might need to use std::tuple to do this, by either subclassing
      // or having calc be a tuple.
      // we would then also need cascade to be a member, rather than directly subclassing.
      // another possibility would be to requires STAGES in the template parameter list
      // and then add calc and stage parameters to set.
      // this might be nice for mixed-use because one might not want all of the Fcs to be the same.
      template<typename T, class F>
      struct df2T : cascade<T, 2, F::stages>
      {
        F calc;
        
        using cascade<T, 2, F::stages>::coeff;
        using cascade<T, 2, F::stages>::state;
        using cascade<T, 2, F::stages>::getCoeffSize;
        using cascade<T, 2, F::stages>::getStateSize;
        
        void set(T omega, T q);

        size_t getStageCount() const { return F::stages; }

        void process(const T* source, T* dest, size_t blockSize);
      };
    };
  }
  
  namespace saturation
  {
    // lovingly borrowed from pichenettes/stmlib
    template<typename T>
    T softlimit(T x) { return x * (27 + x * x) / (27 + 9 * x * x); }

    // lovingly borrowed from pichenettes/stmlib
    template<typename T>
    T softclip(T x) { return x < -3 ? -1 : (x > 3 ? 1 : softlimit(x)); }
  }
  
  struct duration
  {
    // convert bpm to frequency in Hz
    static constexpr double B2F = 1.0 / 60.0;
    static constexpr double F2B = 60;

    // can be used by units an indication for how to treat changes in time values (see: delay)
    enum type : uint8_t
    {
      lerp,
      crossfade,
    };

    double samples; // double so we can express subsample periods.

    // conversions for parameter
    explicit duration(bool b) : samples(b) {}
    explicit duration(int64_t i) : samples(static_cast<double>(i)) {}
    explicit duration(double a) : samples(a) {}
    explicit operator bool() const { return math::abs(samples) >= math::epsilon<double>(); }
    explicit operator int64_t() const { return static_cast<int64_t>(samples); }
    explicit operator double() const { return samples;}

    static duration fromBpm(float bpm, float sampleRate) { return duration(sampleRate/(bpm*B2F)); }
    static duration fromSeconds(float seconds, float sampleRate) { return duration(sampleRate*seconds); }
    float toBpm(float sampleRate) const { return static_cast<float>(F2B*(sampleRate/samples)); }
    float toSeconds(float sampleRate) const { return static_cast<float>(samples/sampleRate); }
  };
  
  // support for tempo detection of a clock signal (i.e. pulse train)
  class clockable
  {
    using period_t = uint32_t;
  protected:
    duration tempo;
    period_t periodMin;
    period_t periodMax;
    period_t ticks;
    float    sampleRate;

    clockable(float sampleRate, period_t samplePeriodMin, period_t samplePeriodMax, float bpm = 60)
    : tempo(duration::fromBpm(bpm, sampleRate)), periodMin(samplePeriodMin), periodMax(samplePeriodMax)
    , ticks(0), sampleRate(sampleRate)
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
    void clock() { tempo.samples = math::constrain(ticks, periodMin, periodMax); ticks = 0; tock(0); }
    void clock(period_t sampleDelay) { tempo.samples = math::constrain<double>(ticks + sampleDelay, periodMin, periodMax); ticks = 0; tock(sampleDelay); }

    float getBpm() const { return tempo.toBpm(sampleRate); }
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
      T evaluate(float phase) const override { return math::sin(math::twoPi<T>() * phase); }
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

  // generates an envelope that begins and ends at zero, with some number of stages leading up to a final stage.
  // each stage of an envelope is defined by a target value, a duration to reach it,
  // and whether the envelope should hold the stage value until it is triggered again.
  template<typename T>
  class envelope : public unitGenerator<T>
  {
  public:
    class stage final : public unitGenerator<T>
    {
      T mTarget;
      T mValue; // current value of the stage
      unit::init<5> init = {
        "stage",
        {
          parameter("target", &mTarget),
          parameter("duration", parameter::type::analog),
          parameter("active", parameter::type::binary),
          parameter("end of stage", parameter::type::binary),
          parameter("value", &mValue)
        }
      };

      T begin; // value the stage started with
      // where we are in the stage
      float time;
      parameter& aw() { return init.params[2]; }
      parameter& ew() { return init.params[3]; }
      
      template<typename E>
      T step();

    public:
      explicit stage(float sampleRate) : unitGenerator<T>(init, sampleRate), mTarget(0), begin(0) { reset(); }

      parameter& target() { return init.params[0]; }
      parameter& duration() { return init.params[1]; }
      const parameter& active() const { return init.params[2]; }
      const parameter& eos() const { return init.params[3]; }
      // current value of the stage
      const parameter& value() const { return init.params[4]; }

      void start(T fromValue) { begin = fromValue; mValue = fromValue; time = -unit::dt(); aw() << true; ew() << false; }
      void reset() { aw() << false; ew() << false; time = -unit::dt(); mValue = 0; }

      template<typename E>
      T generate() { return active() ? step<E>() : mValue; }
      T generate() override { return generate<easing::linear>(); }
    };
  
  private:
    unit::init<1> init = {
      "envelope",
      {
        parameter("eoc", parameter::type::binary)
      }
    };
    array<stage> stages;
    size_t stageIdx;
    stage final;
    parameter& eocw() { return init.params[0]; }

  protected:
    envelope(stage* stages, int stageCount, float sampleRate) : unitGenerator<T>(init, sampleRate)
    , stages(stages, stageCount), stageIdx(0), final(sampleRate) {}
    
  public:
    stage& getStage(size_t idx) { return idx == stages.getSize() ? final : stages[idx]; }
    const stage& getStage(size_t idx) const { return idx == stages.getSize() ? final : stages[idx]; }
    int getStageCount() const { return stages.size() + 1; }

    stage& currentStage() { return getStage(stageIdx); }
    const stage& currentStage() const { return getStage(stageIdx);}
    stage& finalStage() { return final; }
    const stage& finalStage() const { return final; }

    const parameter& value() const { return currentStage().value(); }
    const parameter& eoc() const { return init.params[0]; }

    // make this a parameter we check in generate?
    virtual void trigger();

    template<typename E>
    T generate();
    T generate() override { return generate<easing::linear>(); }

  protected:
    void onSampleRateChanged() override;
    void startStage(int idx, T fromValue) { getStage(idx).start(fromValue); stageIdx = idx; }
    
    // by default, stages advance automatically when their eos goes high.
    // subclasses can override this behavior per stage
    // to enable advancing to the next stage before it is finished,
    // or holding a stage for some period of time.
    virtual bool shouldAdvance(size_t currentStageIdx) { return getStage(currentStageIdx).eos().template read<bool>(); }
  };

  template<typename T>
  class ad : public envelope<T>
  {
    typename envelope<T>::stage attackStage;

  public:
    ad(float attackDuration, float decayDuration, float sampleRate) : envelope<T>(&attackStage, 1, sampleRate), attackStage(sampleRate)
    { attack().target() << T(1); attack().duration() << attackDuration; decay().duration() << decayDuration; }

    typename envelope<T>::stage& attack() { return attackStage; }
    typename envelope<T>::stage& decay() { return envelope<T>::finalStage(); }
    using envelope<T>::eoc;

    using envelope<T>::trigger;
    using envelope<T>::generate;
  };

  template<typename T>
  class asr : public ad<T>
  {
    T triggerThreshold;
    bool gateOn;
    
  public:
    asr(float attackDuration, float decayDuration, float sampleRate, T triggerThreshold = T(0))
    : ad<T>(attackDuration, decayDuration, sampleRate), triggerThreshold(triggerThreshold), gateOn(false) {}

    using ad<T>::attack;
    typename envelope<T>::stage& release() { return envelope<T>::finalStage(); }
    using ad<T>::eoc; 

    void gate(T value);
    void gate(bool on) { gate(on ? T(1) : T(0));}
    void trigger() override { attack().target() << T(1); ad<T>::trigger(); gateOn = false; }
    using envelope<T>::generate;

  protected:
    bool shouldAdvance(size_t currentStageIdx) override { return ad<T>::shouldAdvance(currentStageIdx) && !gateOn; }
  };

  template<typename T>
  class adsr : public envelope<T>
  {
    typename envelope<T>::stage attackStage;
    typename envelope<T>::stage decayStage;
    bool gateOn;

  public:
    adsr(float attackDuration, float decayDuration, float sustainLevel, float releaseDuration, float sampleRate)
    : envelope<T>(&attackStage, 2, sampleRate), attackStage(sampleRate), decayStage(sampleRate), gateOn(false)
    {
      attackStage.duration() << attackDuration;
      attackStage.target() << T(1);
      decayStage.duration() << decayDuration;
      decayStage.target() << sustainLevel;
      envelope<T>::finalStage().duration() << releaseDuration;
    }

    typename envelope<T>::stage& attack() { return attackStage; }
    typename envelope<T>::stage& decay() { return decayStage; }
    parameter& sustain() { return decay().target(); }
    typename envelope<T>::stage& release() { return envelope<T>::finalStage(); }
    using envelope<T>::eoc;

    // should we jump to the release stage if the gate goes off before we start sustaining??
    void gate(bool on)
    {
      if (on && !gateOn)
      {
        envelope<T>::trigger();
      }
      else if (gateOn && !on)
      {
        envelope<T>::startStage(2, envelope<T>::currentStage().value().template read<T>());
      }
      gateOn = on;
    }
    void trigger() override { gateOn = false; envelope<T>::trigger(); }
    using envelope<T>::generate;

  protected:
    bool shouldAdvance(size_t currentStageIdx) override
    {
      return currentStageIdx == 1 ? decayStage.eos() && !gateOn : envelope<T>::shouldAdvance(currentStageIdx);
    }
  };

  template<typename T>
  class slew : public unitProcessor<T>
  {
    T value;
    T epsilon;
    
    unit::init<4> init = {
      "slew",
      {
        parameter("rise", parameter::type::analog),
        parameter("fall", parameter::type::analog),
        parameter("rising", parameter::type::binary),
        parameter("falling", parameter::type::binary)
      }
    };

    parameter& rw() { return init.params[2]; }
    parameter& fw() { return init.params[3]; }

    using unit::dt;
    using unit::getSampleRate;

  public:
    // note: choice of epsilon will depend on the amount of noise in the signal to be slewed.
    // the default value was chosen based on testing with an OWL module's audio input.
    slew(float sampleRate, float riseRate, float fallRate, T initialValue = T(0), T epsilon = math::epsilon<T>()*1000)
    : unitProcessor<T>(init, sampleRate), value(initialValue), epsilon(epsilon)
    { rise() << riseRate; fall() << fallRate; }

    parameter& rise() { return init.params[0]; }
    parameter& fall() { return init.params[1]; }
    const parameter& rising() const { return init.params[2]; }
    const parameter& falling() const { return init.params[3]; }

    T process(const T& v) override;

    using processor<T>::process;
  };

  template<typename T>
  class smoother : public unitProcessor<T>
  {
    T mValue;
    unit::init<2> init = {
      "smoother",
      {
        parameter("degree", parameter::type::analog),
        parameter("value", &mValue)
      }
    };
  public:
    explicit smoother(float smoothingDegree = 0.9f, T initialValue = T(0)) : unitProcessor<T>(init), mValue(initialValue)
    { degree() << smoothingDegree; }

    parameter& degree() { return init.params[0]; }
    const parameter& value() const { return init.params[1]; }

    T process(const T& v) override { return mValue = easing::lerp(v, mValue, math::constrain(*degree(), 0.f, 1.f)); }
    // for block processing
    using unitProcessor<T>::process;
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
    using ring<T>::setWriteIndex;

    // reads behind the write head with sampleDelay (i.e. the ith sample previously written)
    // where a delay of 0 samples will give the most recently written value.
    T read(size_t sampleDelay);

    // reads behind the write head with a fractional sampleDelay and given interpolation
    template<typename I = interpolation::linear<T>>
    T read(float sampleDelay) const;

    // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
    T evaluate(float phase) const override;
  };

  /// delay that will smoothly change time, generating fluctuations in pitch
  template<typename T, typename I = interpolation::linear<T>>
  class delay : public unitProcessor<T>
  {
    duration mTime;
    unit::init<2> init = {
      "delay",
      {
        parameter("time", &mTime),
        parameter("feedback", parameter::type::analog)
      }
    };
    
  protected:
    delayline<T> buffer;
    float mDelayInSamples;

    using unit::dt;
    using unit::getSampleRate;

  public:
    delay(array<T> buffer, float sampleRate, float delayInSeconds = 0, float feedbackAmount = 0)
    : unitProcessor<T>(init, sampleRate), mTime(duration::fromSeconds(delayInSeconds, sampleRate))
    , buffer(buffer.getData(), buffer.getSize())
    {
      mDelayInSamples = mTime.samples;
      feedback() << feedbackAmount;
    }

    /// delay time expressed as vessl::time (i.e. samples), can be set using a double
    parameter& time() { return init.params[0]; }
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    parameter& feedback() { return init.params[1]; }

    T process(const T& in) override;

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<duration::type TimeType = duration::lerp>
    void process(array<T> in, array<T> out);
  };

  template<typename T, typename F>
  class filter : public unitProcessor<T>
  {
    unit::init<2> init = {
      "filter",
      {
        parameter("cutoff", parameter::type::analog),
        parameter("q", parameter::type::analog)
      }
    };
    
    F function;
    float piosr;
    
    using unit::dt;

  public:
    filter(float sampleRate, float cutoffInHz, float kyu = 0) : unitProcessor<T>(init, sampleRate)
    , piosr(math::pi<float>()/sampleRate)
    { cutoff() << cutoffInHz; q() << kyu; }

    parameter& cutoff() { return init.params[0]; }
    parameter& q() { return init.params[1]; }
    
    T process(const T& in) override
    {
      T out;
      function.set(*cutoff()*piosr, *q());
      function.process(&in, &out, 1);
      return out;
    }
    
    void process(array<T> in, array<T> out) override
    {
      function.set(*cutoff()*piosr, *q());
      function.process(in.getData(), out.getData(), in.getSize());
    }

  protected:
    void onSampleRateChanged() override { piosr = math::pi<float>() * dt(); }
  };

  // designed to work with floating point types.
  template<typename T, uint32_t MaxBits>
  class bitcrush : public unitProcessor<T>
  {
    unit::init<3> init = {
      "bit crush",
      {
        parameter("bit rate", parameter::type::analog),
        parameter("bit depth", parameter::type::analog),
        parameter("mangle", parameter::type::binary)
      }
    };
    
    T prevInput;
    T currSample;
    T rateAlpha;

    using unit::dt;

  public:
    bitcrush(float sampleRate, float bitRate, float bitDepth = MaxBits)
      : unitProcessor<T>(init, sampleRate), prevInput(0), currSample(0), rateAlpha(0)
    {
      rate() << bitRate;
      depth() << bitDepth;
    }

    parameter& rate() { return init.params[0]; }
    parameter& depth() { return init.params[1]; }
    parameter& mangle() { return init.params[2]; }

    T process(const T& in) override
    {
      rateAlpha += math::max(1.0f, *rate())*dt();
      if (rateAlpha >= 1)
      {
        rateAlpha -= 1;
        currSample = easing::lerp(prevInput, in, rateAlpha);
      }

      float bd = math::constrain<float>(*depth(), 2.0f, MaxBits);
      float scalar = math::pow(2.f, bd) - 1;
      int val = currSample*scalar;
      if (mangle().template read<bool>())
      {
        val ^= static_cast<int>(prevInput*scalar);
      }
      prevInput = in;
      return static_cast<float>(val) / scalar;
    }

    using unitProcessor<T>::process;
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
  void array<T>::fill(T value)
  {
    writer w (*this);
    while (w)
    {
      w << value;
    }
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

  inline parameter::parameter(const char* name, void* userData): pn(name), pt(type::user)
  {
    pv.u = userData;
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
        bool state = pv.b & 1;
        uint32_t sampleDelay = static_cast<uint32_t>(pv.b >> 1);
        if (sampleDelay == 0) { return T(state); }
        const_cast<parameter*>(this)->pv.b = state | ((sampleDelay - 1) << 1);
        return T(!state);
      }
      case type::digital: return T(pv.i);
      case type::analog: return T(pv.a);
      case type::user: return pv.u ? *static_cast<T*>(pv.u) : T(false);
    }
    return T(false);
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
  
  template<typename T, typename F>
  void filtering::biquad::df2T<T, F>::set(T omega, T q)
  {
    calc(coeff.getData(), omega, q);
  }
  
  template<typename T, typename F>
  void filtering::biquad::df2T<T, F>::process(const T* source, T* dest, size_t blockSize)
  {
    const T* input = source;
    typename array<T>::reader cor(coeff, getCoeffSize());
    for (size_t s = 0; s < getStageCount(); s++)
    {
      T b0 = cor.read();  T b1 = cor.read(); T b2 = cor.read();
      T a1 = cor.read();  T a2 = cor.read();
      T* st = state + 2*s;
      T d1 = st[0]; T d2 = st[1];
      typename array<T>::reader r(input, blockSize);
      typename array<T>::writer w(dest, blockSize);
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

  // @todo ARM specializations for this that utilize the table-based interpolation methods.
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
        static const T RC = static_cast<T>(1) / (math::pi<T> * 200);
        static const T AC = 6.2;
        static T prev = 0;
        float alpha = dt / (dt + RC);
        T white = 2 * (static_cast<T>(random::i32()) / random::I32_MAX) - 1;
        prev = easing::lerp(prev, white, alpha);
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
    if (duration() > math::epsilon<float>())
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
    return easing::lerp(mFrom, mTo, lt);
  }

  template<typename T>
  template<typename E>
  T envelope<T>::stage::step()
  {
    float dt = unit::dt();
    float s = unit::dt() / math::max(*duration(), dt);
    time += s;
    float t = math::constrain(time, 0.f, 1.f);
    mValue = easing::interp<E, T>(begin, mTarget, t);
    if (time >= 1)
    {
      aw() << false;
      ew() << true;
    }
    return mValue;
  }

  template<typename T>
  void envelope<T>::trigger()
  {
    for (stage& stage : stages)
    {
      stage.reset();
    }
    final.reset();
    eocw() << false;
    stageIdx = 0;
    stages[0].start(0);
  }

  template<typename T>
  template<typename E>
  T envelope<T>::generate()
  {
    T value = currentStage().template generate<E>();
    if (stageIdx == stages.getSize() && final.eos())
    {
      eocw() << true;
    }
    else if (shouldAdvance(stageIdx))
    {
      getStage(++stageIdx).start(value);
    }
    return value;
  }

  template<typename T>
  void envelope<T>::onSampleRateChanged()
  {
    for (stage& stage : stages)
    {
      stage.setSampleRate(unit::getSampleRate());
    }
    final.setSampleRate(unit::getSampleRate());
  }

  template<typename T>
  void asr<T>::gate(T value)
  {
    bool valueOn = value > triggerThreshold;
    T attackTarget = attack().target().template read<T>();
    if (!gateOn && valueOn)
    {
      ad<T>::trigger();
      gateOn = true;
      // reset the attackTarget because where the gate value ends up might be less than the previous attackTarget
      attackTarget = 0;
    }
    else if (gateOn && !valueOn)
    {
      gateOn = false;
      // if we're mid-attack, start the release
      if (attack().active())
      {
        envelope<T>::startStage(1, attack().value().template read<T>());
      }
    }
    attack().target() << math::max(value, attackTarget);
  }

  template<typename T>
  T slew<T>::process(const T& v)
  {
    bool isRise = v > value+epsilon;
    bool isFall = v < value-epsilon;
    float rate = isRise ? *rise() : isFall ? -1.0f*fall() : (v-value)*getSampleRate();
    value += rate*dt();
    rw() << isRise;
    fw() << isFall;
    return value;
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
    assert(sampleDelay < getSize());
    sampleDelay = getSize() - 1 - sampleDelay;
    size_t idx = getWriteIndex() + sampleDelay;
    return getData()[idx % getSize()];
  }

  template<typename T>
  template<typename I>
  T delayline<T>::read(float sampleDelay) const
  {
    assert(sampleDelay >= 0 && sampleDelay <= getSize()-1);
    float fSize = static_cast<float>(getSize());
    sampleDelay = fSize - 1 - sampleDelay;
    float fidx = static_cast<float>(getWriteIndex()) + sampleDelay;
    float idx;
    float f = math::mod(fidx, &idx);
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
    mDelayInSamples = easing::lerp<float>(mDelayInSamples, mTime.samples, dt() * 10);
    // delay time in samples
    float dts = math::constrain(mDelayInSamples, 0.f, static_cast<float>(buffer.getSize()-1));
    float s = buffer.template read<I>(dts);
    float fbk = math::constrain(feedback() >> fbk, -1.f, 1.f);
    buffer.write(in + s * fbk);
    return s;
  }

  template<typename T, typename I>
  template<duration::type TimeType>
  void delay<T, I>::process(array<T> in, array<T> out)
  {
    if (TimeType == duration::crossfade)
    {
      float fade = 0;
      float fadeInc = 1.0f / in.getSize();
      // smooth time parameter to prevent crunchiness when it is noisy or changes by large amounts
      float targetSampleDelay = easing::lerp<float>(mDelayInSamples, mTime.samples, dt() * 10);
      // delay time in samples
      float fts = math::constrain(mDelayInSamples, 0.f, static_cast<float>(buffer.getSize()-1));
      float tts = math::constrain(targetSampleDelay, 0.f, static_cast<float>(buffer.getSize()-1));
      float fbk = math::constrain(*feedback(), -1.f, 1.f);

      typename array<T>::reader r = in.getReader();
      typename array<T>::writer w = out.getWriter();
      while (r)
      {
        float wet = (1.0f - fade) * buffer.template read<I>(fts) + fade * buffer.template read<I>(tts);
        buffer.write(r.read() + wet*fbk);
        w << wet;
        fade += fadeInc;
      }
        
      mDelayInSamples = targetSampleDelay;
    }
    else
    {
      processor<T>::process(in, out);
    }
  }

  template<>
  inline float math::wrap01(float v) { float i; return modf(v, &i); }

#ifdef ARM_CORTEX
  template<>
  void array<float>::fill(float value)
  {
    arm_fill_f32(value, data, size);  
  }

  template<>
  array<float> array<float>::operator<<(array copyFrom)
  {
    VASSERT(this->getSize() >= copyFrom.getSize(), "Not enough room in this array for contents of copyFrom");
    arm_copy_f32(copyFrom.getData(), this->getData(), copyFrom.getSize());
    return *this;
  }
  
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
  
  template<>
  inline float math::sin(float x) { return arm_sin_f32(x); }

  template<>
  inline float math::sqrt(float x)
  {
    float out;
    if (ARM_MATH_SUCCESS == arm_sqrt_f32(x, &out))
    {
      return out;
    }
    return 0;
  }
  
  template<typename F>
  struct filtering::biquad::df2T<float, F> : cascade<float, 2, F::stages>
  {
    F calc;
    arm_biquad_cascade_df2T_instance_f32 inst;
          
    using cascade<float, 2, F::stages>::coeff;
    using cascade<float, 2, F::stages>::state;
    using cascade<float, 2, F::stages>::getCoeffSize;
    using cascade<float, 2, F::stages>::getStateSize;

    df2T() { arm_biquad_cascade_df2T_init_f32(&inst, getStageCount(), coeff.getData(), state.getData()); }
          
    void set(float omega, float q)
    {
      calc(coeff.getData(), omega, q);
    }

    size_t getStageCount() const { return F::stages; }

    void process(const float* source, float* dest, size_t blockSize)
    {
      arm_biquad_cascade_df2T_f32(&inst, source, dest, blockSize);
    }
  };
#endif
}
