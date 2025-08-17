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
#include <cmath>
#include <cfloat>
#include <cassert>
#include <cstring> // for memcpy

// When built for ARM Cortex-M processor series,
// we provide template specializations that use the optimized CMSIS library:
// http://www.keil.com/pack/doc/CMSIS/General/html/index.html
#ifdef ARM_CORTEX
#include "arm_math.h" 
#endif //ARM_CORTEX

#ifndef PI
#define PI 3.14159265358979323846
#endif

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
  class processor
  {
  public:
    processor() = default;
    virtual ~processor() = default;
    processor(const processor&) = delete;
    processor(const processor&&) = delete;
    processor& operator=(const processor&) = delete;
    processor& operator=(processor&&) = delete;
    
    virtual T process(const T& in) = 0;  // NOLINT(portability-template-virtual-member-function)
    
    void process(source<T>& in, sink<T>& out)
    {
      while (in && out)
      {
        out.write(process(in.read()));
      }
    }
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

    class reader : public source<T>
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

    class writer : public sink<T>
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

    array operator<<(array copyFrom)
    {
      writer w(*this);
      w << reader(copyFrom);
      return *this;
    }

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
  T* begin(array<T>& arr)
  {
    return arr.begin();
  }

  template<typename T>
  T* end(array<T>& arr)
  {
    return arr.end();
  }

  template<typename T>
  typename array<T>::writer array<T>::writer::operator<<(const reader& r) {
    size_t rsz = r.available();
    size_t wsz = available();
    assert(wsz >= rsz);
    const T* rh = *r;
    memcpy(static_cast<void*>(head), static_cast<const void*>(rh), rsz*sizeof(T));
    head += rsz;
    return writer(head, wsz - rsz);
  }
  
  template<typename T>
  array<T> array<T>::add(array other, array dest)
  {
    assert(size == other.size && size <= dest.size);
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
    assert(size <= dest.size);
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
  class ring : array<T>
  {
    T* head;
  public:
    ring(T* inData, size_t inSize) : array<T>(inData, inSize), head(inData + inSize - 1)
    {
    }

    // expose direct access to underlying array data
    using array<T>::getData;
    using array<T>::getSize;

    size_t getWriteIndex() const
    {
      return head - array<T>::data;
    }
    
    void write(const T& v)
    {
      *head++ = v;
      if (head == array<T>::end())
      {
        head = array<T>::begin();
      }
    }

    ring operator<<(typename array<T>::reader r)
    {
      assert(r.available() < getSize());
      while (r)
      {
        write(r.read());
      }
      return *this;
    }
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

    parameter(const char * name, type type) : pn(name), pt(type)
    {
      switch (pt)
      {
        case type::binary: pv.b = false; break;
        case type::digital: pv.i = 0; break;
        case type::analog: pv.a = 0; break;
        case type::user: pv.u = nullptr; break;
          
        default: assert(false && "attempted to initialize a parameter with an unknown type"); break;  // NOLINT(clang-diagnostic-covered-switch-default)
      }
    }
    parameter(const char * name, bool value) : pn(name), pt(type::binary) { pv.b = value; }
    parameter(const char * name, long value) : pn(name), pt(type::digital) { pv.i = value; }
    parameter(const char * name, double value) : pn(name), pt(type::analog) { pv.a = value; }
    parameter(const char * name, void* value) : pn(name), pt(type::user) { pv.u = value; } 

    const char* getName() const { return pn; }
    type getType() const { return pt; }

    // for the sane
    template<typename T>
    T read() const
    {
      // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
      switch (pt)
      {
        case type::binary: return static_cast<T>(pv.b);
        case type::digital: return static_cast<T>(pv.i);
        case type::analog: return static_cast<T>(pv.a);
        case type::user: return pv.u ? *static_cast<T*>(pv.u) : T();
      }
      return T();
    }

    template<typename T>
    parameter& write(const T& value)
    {
      // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
      switch (pt)
      {
        case type::binary: pv.b = static_cast<bool>(value); break;
        case type::digital: pv.i = static_cast<long>(value); break;
        case type::analog: pv.a = static_cast<double>(value); break;
        case type::user: if (pv.u) *static_cast<T*>(pv.u) = value; break;
      }
      return *this;
    }

    // overloading dereference with float conversion because it will be used so often
    float operator*() const { return read<float>(); }

    template<typename T>
    T operator!() const
    {
      // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
      switch (pt)
      {
        case type::binary: return !pv.b;
        case type::digital: return !pv.i;
        case type::analog: return !pv.a;
        case type::user: return !(pv.u ? *pv.u : T());
      }
      return !T(); 
    }

    template<typename T>
    T operator-() const
    {
      // ReSharper disable once CppDefaultCaseNotHandledInSwitchStatement
      switch (pt)
      {
        case type::binary: return !pv.b;
        case type::digital: return -pv.i;
        case type::analog: return -pv.a;
        case type::user: return -(pv.u ? *pv.u : T());
      }
      return -T(); 
    }

  private:
    const char* pn;
    union
    {
      bool b;
      long i;
      double a;
      void* u;
    } pv;
    type pt;
  };

  template<typename T> parameter& operator<<(parameter& p, const T& v) { return p.write(v); }
  template<typename T> T& operator>>(parameter& p, T& v) { v = p.read<T>(); return v; }
  template<typename T> bool operator>(const parameter& p, const T& v) { return p.read<T>() > v; }
  template<typename T> bool operator<(const parameter& p, const T& v) { return p.read<T>() < v; }
  template<typename T> bool operator==(const parameter& p, const T& v) { return p.read<T>() == v; }
  template<typename T> bool operator!=(const parameter& p, const T& v) { return p.read<T>() != v; }
  template<typename T> T operator+(const parameter& p, const T& v) { return p.read<T>() + v; }
  template<typename T> T operator+(const T& v, const parameter& p) { return p.read<T>() + v; }
  template<typename T> T operator*(const parameter& p, const T& v) { return p.read<T>() * v; }
  template<typename T> T operator*(const T& v, const parameter& p) { return p.read<T>() * v; }
  template<typename T> T operator/(const parameter& p, const T& v) { return p.read<T>() / v; }
  template<typename T> T operator/(const T& v, const parameter& p) { return v / p.read<T>(); }
  
  class unit : array<parameter>
  {
  protected:
    template<size_t N>
    struct init
    {
      const char* name;
      parameter params[N];
    };
    
    template<size_t N>
    explicit unit(init<N>& init, float sampleRate = 1) : array(init.params, N), n(init.name) { setSampleRate(sampleRate); }
    unit(const char* name, parameter* params, int paramsCount, float sampleRate = 1) : array(params, paramsCount), n(name) { setSampleRate(sampleRate); }

  public:
    virtual ~unit() = default;
    unit(const unit&) = default;
    unit(unit&&) = default;
    unit& operator=(const unit&) = default;
    unit& operator=(unit&&) = default;
    
    const char* name() const { return n; }
    
    float getSampleRate() const { return sampleRate; }
    
    void setSampleRate(float sr)
    {
      sampleRate = sr;
      deltaTime = 1.0f/sr;
      onSampleRateChanged();
    }

    parameter& getParameter(size_t index) { return this->operator[](index); }
    const parameter& getParameter(size_t index) const { return this->operator[](index); }
    size_t getParameterCount() const { return getSize(); }

    using array::begin;
    using array::end;

  protected:
    virtual void onSampleRateChanged() {}
    float dt() const { return deltaTime; }
    
  private:
    const char* n;
    float sampleRate;
    float deltaTime;
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
      T operator()(const T* buffer, double fracIdx)
      {
        return buffer[static_cast<size_t>(round(fracIdx))];
      }
    };

    template<typename T>
    struct linear
    {
      T operator()(const T* buffer, float fracIdx)
      {
        float idx;
        float frac = modf(fracIdx, &idx);
        size_t x0 = static_cast<size_t>(idx);
        return buffer[x0] + (buffer[x0 + 1] - buffer[x0])*frac;
      }
    };

    template<typename T>
    struct cubic
    {
      T operator()(const T* buffer, float fracIdx)
      {
        static const T DIV6 = static_cast<T>(1. / 6.);
        static const T DIV2 = static_cast<T>(0.5);

        float idx;
        float f = modf(fracIdx, &idx);
        float fm1 = f - 1.f;
        float fm2 = f - 2.f;
        float fp1 = f + 1.f;
        size_t x0 = idx;
        return -f * fm1*fm2*DIV6 * buffer[x0 - 1] + fp1 *fm1*fm2*DIV2 * buffer[x0] - fp1 * f*fm2*DIV2 * buffer[x0 + 1] + fp1 * f*fm1*DIV6 * buffer[x0 + 2];
      }
    };
  };
  
  template<typename T, typename I = interpolation::linear<T>>
  T sample(const T* buffer, float fracIdx)
  {
    assert(fracIdx >= 0 && "fracIdx argument to sample must be non-negative!");
    static I interpolator;
    return interpolator(buffer, fracIdx);
  }

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

  template<>
  inline float wrap01(float v) { float i; return modf(v, &i); }

#ifndef clamp
  template<typename T>
  T clamp(T val, T low, T high)
  {
    return val < low ? low : val > high ? high : val;
  }
#endif

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
      T evaluate(float phase) const override { return sin(2*PI*phase); }
    };
    
    template<>
    inline float sine<float>::evaluate(float phase) const { return sin(static_cast<float>(2 * PI * phase)); }
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
    explicit wavetable(source<T>& source): waveform<T>()
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

    explicit wavetable(const waveform<T>& waveform): waveform<T>()
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

    // ReSharper disable once CppMemberFunctionMayBeStatic
    size_t size() const
    {
      return N;
    }

    T get(const size_t i) const
    {
      return buffer[i+1];
    }

    void set(const size_t i, T val)
    {
      buffer[i + 1] = val;
      // also update wrap-around values
      switch (i)
      {
        case 0: buffer[N+1] = buffer[1]; break;
        case 1: buffer[N+2] = buffer[2]; break;
        case N - 1: buffer[0] = buffer[N]; break;
        default: break;
      }
    }
    
    // implement waveform:
    T evaluate(float phase) const override
    {
      float idx = phase*N;
      while (idx > N) idx -= N;
      while (idx < 0) idx += N;
      return sample<T, I>(buffer, idx);
    }
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
    T read(size_t sampleDelay)
    {
      assert(sampleDelay < getSize()-1);
      size_t idx = getWriteIndex() + 1 + sampleDelay;
      return getData[idx%getSize()];
    }

    // reads behind the write head with a fractional sampleDelay and given interpolation
    template<typename I = interpolation::linear<T>>
    T read(float sampleDelay) const
    {
      float fSize = static_cast<float>(getSize());
      sampleDelay = fSize - sampleDelay;
      assert(sampleDelay >= 0 && sampleDelay < fSize - 1);
      sampleDelay = static_cast<float>(getWriteIndex() + 1) + sampleDelay;
      float idx;
      float f = modf(sampleDelay, &idx);
      size_t x0 = static_cast<size_t>(idx) % getSize();
      size_t x1 = (x0 + 1) % getSize();
      size_t x2 = (x0 + 2) % getSize();
      const T* data = getData();
      T s[3] = { data[x0], data[x1], data[x2] };
      return sample<T, I>(s, f);
    }
    
    // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
    T evaluate(float phase) const override
    {
      float fSize = static_cast<float>(getSize());
      phase = wrap(phase, -1.f, 1.f);
      float sampleDelay = phase > 0 ? (1.0f - phase)*fSize : -phase*fSize;
      return read(sampleDelay);
    }
  };

  // @todo easings namespace similar to interpolation namespace
  template<typename T>
  T lerp(T v1, T v2, float a)
  {
    return v1 + (v2 - v1)*a;
  }

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
      "noise", { parameter("tint", parameter::type::digital), parameter("rate", parameter::type::analog) }
    };

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
    
    T generate() override
    {
      step += dt() * rate().read();
      float alpha = step / dt();
      if (alpha >= 1)
      {
        noiseTint tnt = static_cast<noiseTint>(tint());
        nz[0] = nz[1];
        nz[1] = next(tnt, dt());
        alpha = wrap01(alpha);
        step = alpha * dt();
      }
      return lerp(nz[0], nz[1], alpha);
    }

  private:
    float step;
    T nz[2];
  };
  
  template<typename T>
  T noiseGenerator<T>::next(noiseTint tint, double dt)
  {
    switch (tint)
    {
      case noiseTint::white:
      {
        return 2 * (static_cast<T>(rand()) / RAND_MAX) - 1;
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
        T white = 2 * (static_cast<T>(rand()) / RAND_MAX) - 1;
        prev = lerp<T>(prev, white, alpha);
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
        static int whiteValues[6] = { rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6) };

        int lastKey = key;
        T sum = 0;
        key = key == MAX_KEY ? 0 : ++key;

        int diff = lastKey ^ key;
        for (int i = 0; i < 6; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            whiteValues[i] = rand() % (RANGE / 6);
          }
          sum += whiteValues[i];
        }
        maxSum = std::max<T>(sum, maxSum);
        T n = 2 * (sum / maxSum) - 1;
        assert(!isnan(n) && "pink noise generated nan");
        return n;
      }
    }

    return 0;
  }

  
  // unit that generates a linear ramp from one value to another over a duration of seconds
  // @todo implement easings above and add that as a template parameter
  template<typename T>
  class ramp final : public unitGenerator<T>
  {
    T mFrom, mTo;
    unit::init<4> init {
      "ramp",
      { parameter("from", &mFrom), parameter("to", &mTo), parameter("duration", parameter::type::analog), parameter("eor", parameter::type::binary)}
    };

    using unit::dt;
    
  public:
    explicit ramp(float sampleRate, float durationInSeconds = 0, T fromValue =T(0), T toValue = T(0)) : unitGenerator<T>(init, sampleRate)
    , t(0)
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
  
    void trigger()
    {
      if (duration() > FLT_EPSILON)
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
  
    T generate() override
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
      return lerp<T>(mFrom, mTo, lt);
    }
  
  private:
    parameter& eorw() { return init.params[3]; }
    float t;
  };
  //
  // template<typename T, int I, int O>
  // class mixer final : public unit<T>
  // {
  // public:
  //   mixer(T masterVolume)
  //     : audio(this)
  //     , volProxy(volCtrl)
  //     , masterCtrl(this, masterVolume)
  //   {
  //     // default all volume controls to 1
  //     // and the mix matrix to "unity",
  //     // meaning that each input is sent only to the corresponding output at full volume.
  //     for (int i = 0; i < I; ++i)
  //     {
  //       vol[i] = 1;
  //       for (int o = 0; o < O; ++o)
  //       {
  //         matrix[i][o] = T(i == o);
  //       }
  //     }
  //   }
  //
  //   /// audio inputs
  //   array<input>& in = audio.inputs();
  //   // audio outputs
  //   array<output>& out = audio.outputs();
  //
  //   /// mix matrix - indicates how much of each input signal to send to each output
  //   T matrix[I][O];
  //
  //   /// volume modulation for each input
  //   array<input>& vol = volProxy;
  //
  //   /// master volume applied to each output after accumulating inputs
  //   input& master = masterCtrl;
  //
  //   array<input>& inputs() override { return in; }
  //   array<output>& outputs() override { return out; }
  //
  //   void tick(T deltaTime) override
  //   {
  //     T mv = masterCtrl;
  //     for (int o = 0; o < O; ++o)
  //     {
  //       // initialize
  //       T v = 0;
  //       // sum all inputs to this output based on matrix and volume settings
  //       for (int i = 0; i < I; ++i)
  //       {
  //         v += audio.in[i] * matrix[i][o] * volCtrl[i] * mv;
  //       }
  //       audio.out[o] = v;
  //     }
  //   }
  //
  // private:
  //   // audio in/out
  //   io<T, I, O> audio;
  //   value<T> volCtrl[I];
  //   proxyarray<T, I, input> volProxy;
  //   value<T> masterCtrl;
  // };
  //
  // // lightweight unit to transform an input value to an output via a user-provided delegate
  // template<typename T, int I = 1>
  // class function : public unit<T>
  // {
  // public:
  //   using delegateType = std::function<T(T, T)>;
  //
  //   function(delegateType d)
  //     : io(this)
  //     , delegate(d)
  //   {
  //   }
  //
  //   delegateType delegate;
  //
  //   input& in = io.in[0];
  //   output& out = io.out[0];
  //
  //   array<input>& inputs() override { return io.inputs(); }
  //   array<output>& outputs() override { return io.outputs(); }
  //
  //   void tick(T deltaTime) override
  //   {
  //     io.out[0] = delegate(io.in[0], deltaTime);
  //   }
  //
  // protected:
  //   io<T, I, 1> io;
  // };
  //
  // // multiplies 'in' by 'factor' and puts the result in 'out'
  // template<typename T>
  // class multiplier final : public function<T, 2>
  // {
  // public:
  //    multiplier(T amount = 0) 
  //     : function([this](T val, T dt){ return val * io.in[1]; })
  //   {
  //     io.in[1] = amount;
  //   }
  //
  //   unit::input& factor = io.in[1];
  // };
  //
  // // sums all inputs to a single output
  // template<typename T, int I>
  // class summer final : public unit<T>
  // {
  // public:
  //   summer() : io(this) { }
  //
  //   array<input>& in = io.inputs();
  //   output& out = io.out[0];
  //
  //   array<input>& inputs() override { return in; }
  //   array<output>& outputs() override { return io.outputs(); }
  //
  //   void tick(T deltaTime) override
  //   {
  //     T sum = 0;
  //     for (int i = 0, sz = in.size(); i < sz; ++i)
  //     {
  //       sum += io.in[i];
  //     }
  //     io.out[0] = sum;
  //   }
  //
  // private:
  //   io<T, I, 1> io;
  // };

  template<typename T>
  class oscil final : public unitGenerator<T>
  {
    unit::init<4> init = {
      "oscil", {
        parameter("frequency", parameter::type::analog), parameter("fm (lin)", parameter::type::analog),
        parameter("fm (v/oct)", parameter::type::analog), parameter("phase mod", parameter::type::analog)
      }
    };

    waveform<T>* wave;
    float phase;
    
    using unit::dt;
    
  public:
    explicit oscil(float sampleRate, waveform<T>& wave, double freqInHz = 440) : unitGenerator<T>(init, sampleRate)
    , wave(&wave), phase(0) 
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
    
    T generate() override
    {
      T val = wave->evaluate(phase + pm());
      phase += (fHz() * exp2(*fmExp()) + fmLin())*dt();
      phase = wrap01(phase);
      return val;
    }

    void reset()
    {
      phase = 0;
    }
  };
  
//
//   template<typename T>
//   class waveshaper final : public unit<T>
//   {
//   public:
//     using waveform = waveform<T>;
//
//     waveshaper(const waveform& waveShape, bool wrapInput = true)
//       : shape(waveShape)
//       , wrap(wrapInput)
//       , io(this)
//     {
//
//     }
//
//     // input in [-1,1] range that will be used sample the shape waveform
//     input& in = io.in[0];
//     // sampled shape value based on current value of in
//     output& out = io.out[0];
//
//     const waveform& shape;
//     bool wrap = true;
//
//     array<input>& inputs() override { return io.inputs(); }
//     array<output>& outputs() override { return io.outputs(); }
//
//     void tick(T deltaTime) override
//     {
//       assert(!isnan(io.in[0]) && "waveshaper input is nan!");
//       T phase = wrap ? wrapPhase<T>(io.in[0] * 0.5 + 0.5) : clamp<T>(io.in[0] * 0.5 + 0.5, 0, 1);
//       io.out[0] = shape.value(phase);
//     }
//
//   private:
//     io<T, 1, 1> io;
//   };

   template<typename T, typename I = interpolation::linear<T>>
   class delay : public unitProcessor<T>
   {
     float mTime;
     unit::init<2> init = {
       "delay", { parameter("time", parameter::type::analog), parameter("feedback", parameter::type::analog)}
     };
     delayline<T> buffer;

     using unit::dt;
   
   public:
     delay(array<T> buffer, float sampleRate, float delayInSeconds = 0, float feedbackAmount = 0) : unitProcessor<T>(init, sampleRate)
     , mTime(delayInSeconds), buffer(buffer.getData(), buffer.getSize())
     {
       time() << delayInSeconds;
       feedback() << feedbackAmount;
     }
     
     /// amount of delay time in seconds
     parameter& time() { return init.params[0]; }
     /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
     parameter& feedback() { return init.params[1]; }

     T process(const T& in) override
     {
       // smooth time parameter to prevent crunchiness when it is noisy or changes by large amounts
       mTime = lerp(mTime, *time(), dt()*10);
       // delay time in samples
       float dts = clamp(mTime * unit::getSampleRate(), 0, (float)buffer.getSize()-1);
       float fbk = clamp(feedback() >> fbk, -1.f, 1.f);
       T s = buffer.template read<I>(dts);
       buffer.write(in + s*fbk);
       return s;
     }

     using processor<T>::process;
   };

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

// need to break out the rotation matrix calculation into a separate class
// to allow for easier template specialization.
// but, since rotation matrices beyond 3x3 start to get real math-y,
// I'm not even sure this is something that should be included in this class.
// it might make more sense to have a separate class that is specifically about
// positioning inputs in 2d or 3d "space" parametrically,
// at let this class be a generic tool for routing inputs to outputs with scaling.
//void calcRotationMatrix()
//{
//  T rad = rotCtrl.in[0] * PI;
//  T sin = std::sin(rad);
//  T cos = std::cos(rad);
//  rotation[0][0] = cos; rotation[0][1] = -sin;
//  rotation[1][0] = sin; rotation[1][1] = cos;
//}
//
//for (int i = 0; i < I; ++i)
//{
//  T* mix = matrix[i];
//  for (int o = 0; o < O; ++o)
//  {
//    // for every in -> out mix coefficient m
//    // we need to compute m-prime,
//    // by performing a matrix multiplication across
//    // the rotation matrix's o-th row.
//    float mp = 0;
//    int ri = o;
//    for (int ro = 0; ro < O; ++ro)
//    {
//      mp += mix[ro] * rotation[ri][ro];
//    }
//
//    // we then scale the rotated mix coefficient by overall volume and apply to i-th input
//    io.out[o] = v * mp*io.in[i];
//  }
//}
