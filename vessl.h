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
#include <array>
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

// Note: In all classes using typename T, it is assumed to be POD and to have support for all arithmetic operators
namespace vessl
{
  using char_t = char;
  using size_t = std::size_t;
  using binary_t = bool;
  using digital_t = int64_t;
  using analog_t = float;
  
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

    virtual binary_t isEmpty() const = 0;
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

    virtual binary_t isFull() const = 0;
    virtual void write(const T& value) = 0;

    explicit operator binary_t() const { return !isFull(); }
  };

  template<typename T>
  class generator : public source<T>
  {
  public:
    generator() = default;
    virtual T generate() = 0;

    // by default, we assume an endless source of data
    binary_t isEmpty() const override { return false; }
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
    binary_t isEmpty() const { return size == 0; }
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
      binary_t isEmpty() const override { return head == end; }
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

      binary_t isFull() const override { return head == end; }
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
    array add(array other, array dest) const;
    // returns this
    array add(array other) { return add(other, *this); }
    array subtract(array other, array dest) const;
    array subtract(array other) { return subtract(other, *this); }
    // returns dest
    array scale(analog_t value, array dest) const;
    // returns this
    array scale(analog_t value) { return scale(value, *this); }
    // returns dest
    array multiply(array other, array dest) const;
    // returns this
    array multiply(array other) { return multiply(other, *this); }
  };

  template<typename T>
  T* begin(array<T>& arr) { return arr.begin(); }

  template<typename T>
  T* end(array<T>& arr) { return arr.end(); }
  
  namespace frame
  {
    template<typename T, size_t N>
    struct channels;
    
    namespace mono
    {
      typedef channels<analog_t, 1> analog_t;
      typedef channels<digital_t, 1> digital_t;
    }

    namespace stereo
    {
      typedef channels<analog_t, 2> analog_t;
      typedef channels<digital_t, 2> digital_t;
    }
    
    template<typename T, size_t N>
    struct channels : array<T>
    {
      T samples[N];
      channels() : array<T>(samples, N) { array<T>::fill(0); }
      explicit channels(T m) : array<T>(samples, N) { array<T>::fill(m); }
      channels(channels const& other) : array<T>(samples, N) { array<T>::copy(other); }
      // @todo move and assignment operators
      
      channels<T, 1> toMono() const
      {
        T sum = 0;
        for (size_t c = 0; c < N; ++c)
        {
          sum += samples[c];
        }
        return channels<T, 1>(sum / N);
      }
    };

    template<typename T>
    struct channels<T, 1> : array<T>
    {
      T samples[1];
      channels() : array<T>(samples, 1) { samples[0] = 0; }
      explicit channels(T m) : array<T>(samples, 1) { samples[0] = m; }
      channels(const channels& other) : array<T>(samples, 1) { samples[0] = other.samples[0]; }
      channels(channels&& other) noexcept : array<T>(samples, 1) { samples[0] = std::move(other.samples[0]); }
      virtual ~channels() = default;
      
      channels& operator=(const channels& other)  // NOLINT(bugprone-unhandled-self-assignment)
      {
        if (this == &other)
        {
          return *this;
        }
        
        samples[0] = other.samples[0];
        return *this;
      }
      
      channels& operator=(channels&& other) noexcept
      {
        if (this == &other)
        {
          return *this;
        }
        
        samples[0] = std::move(other.samples[0]);
        return *this;
      }
      
      channels toMono() const;

      T& value() { return samples[0]; }
      const T& value() const { return samples[0]; }
    };

    template<typename T>
    struct channels<T, 2> : array<T>
    {
      T samples[2];

      channels() : array<T>(samples, 2) { samples[0] = 0; samples[1] = 0; }
      explicit channels(T m) : array<T>(samples, 2) { samples[0] = m; samples[1] = m; }
      channels(T left, T right) : array<T>(samples, 2) { samples[0] = left, samples[1] = right; }
      channels(const channels& other) : array<T>(samples, 2) { samples[0] = other.samples[0]; samples[1] = other.samples[1]; }
      channels(channels&& other) noexcept : array<T>(samples, 2) { samples[0] = std::move(other.samples[0]); samples[1] = std::move(other.samples[1]); }
      virtual ~channels() = default;
      
      channels& operator=(const channels& other)  // NOLINT(bugprone-unhandled-self-assignment)
      {
        if (this == &other)
        {
          return *this;
        }
        
        samples[0] = other.samples[0];
        samples[1] = other.samples[1];
        return *this;
      }
      
      channels& operator=(channels&& other) noexcept
      {
        if (this == &other)
        {
          return *this;
        }
        
        samples[0] = std::move(other.samples[0]);
        samples[1] = std::move(other.samples[1]);
        return *this;
      }
      
      channels<T, 1> toMono() const { return channels<T, 1>((samples[0] + samples[1]) * 0.5f); }
    
      T& left() { return samples[0]; }
      const T& left() const { return samples[0]; }
      T& right() { return samples[1]; }
      const T& right() const { return samples[1]; }
    };

    template<typename T, size_t N>
    constexpr channels<T, N> operator+(channels<T,N> lhs, const channels<T,N>& rhs)
    {
      channels<T, N> result;
      lhs.add(rhs, result);
      return result;
    }

    template<typename T, size_t N>
    constexpr channels<T, N> operator-(channels<T,N> lhs, const channels<T,N>& rhs)
    {
      channels<T, N> result;
      lhs.subtract(rhs, result);
      return result;
    }

    template<typename T, size_t N>
    constexpr channels<T, N> operator*(channels<T, N> lhs, const channels<T,N>& rhs)
    {
      channels<T, N> result;
      lhs.multiply(rhs, result);
      return result;
    }

    template<typename T, size_t N>
    constexpr channels<T, N> operator*(channels<T, N> lhs, const analog_t& rhs)
    {
      channels<T, N> result;
      lhs.scale(rhs, result);
      return result;
    }

    template<typename T, size_t N>
    constexpr channels<T, N> operator*(analog_t lhs, const channels<T, N>& rhs)
    {
      channels<T, N> result;
      rhs.scale(lhs, result);
      return result;
    }

    template<typename T, size_t N>
    constexpr channels<T, N> operator^(channels<T, N> lhs, const channels<T,N>& rhs)
    {
      channels<T, N> result;
      for (size_t i = 0; i < N; ++i)
      {
        result[i] = static_cast<digital_t>(lhs[i]) ^ static_cast<digital_t>(rhs[i]);
      }
      return result;
    }
  }

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
      binary = 0, // on/off (binary_t)
      digital = 1, // integral values (digital_t)
      analog = 2, // floating point values (analog_t)
      // space for more built-ins

      // user provided type, stored as a void*, must be convertible to all other parameter types.
      // even if the conversion is meaningless.
      user = UINT8_MAX 
    };

    parameter(const char_t* name, type type);
    parameter(const char_t* name, void* userData);

    const char_t* getName() const { return pn; }
    type getType() const { return pt; }

    // for the sane, static cast from our stored type to T
    template<typename T>
    T read() const;

    binary_t readBinary() const { return read<binary_t>(); }
    explicit operator binary_t() const { return read<binary_t>(); }
    
    // read the actual state of a binary parameter, including sample delay
    // @todo rename to readBinaryRaw
    binary_t read(uint32_t* sampleDelay) const;

    digital_t readDigital() const { return read<digital_t>(); }
    analog_t readAnalog() const { return read<analog_t>(); }

    // overloading dereference with analog conversion because it will be used so often
    analog_t operator*() const { return read<analog_t>(); }

    template<typename T>
    parameter& write(const T& value);

    // set the state of a binary parameter with a sample delay
    parameter& write(binary_t gate, uint32_t sampleDelay);

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

  private:
    const char_t* pn;

    union
    {
      size_t  b;
      digital_t i;
      analog_t  a;
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
    analog_t sampleRate;
    analog_t deltaTime;
    
  protected:
    template<size_t N>
    struct init
    {
      const char_t* name;
      parameter params[N];
    };

    template<size_t N>
    explicit unit(init<N>& init, analog_t sampleRate = 1) : array(init.params, N), unitName(init.name) { setSampleRate(sampleRate); }
    unit(const char_t* name, parameter* params, size_t paramsCount, analog_t sampleRate = 1) : array(params, paramsCount), unitName(name) { setSampleRate(sampleRate); }

  public:
    virtual ~unit() = default;
    unit(const unit&) = default;
    unit(unit&&) = default;
    unit& operator=(const unit&) = default;
    unit& operator=(unit&&) = default;

    const char_t* name() const { return unitName; }

    analog_t getSampleRate() const { return sampleRate; }

    void setSampleRate(analog_t sr);

    parameter& getParameter(size_t index) { return this->operator[](index); }
    const parameter& getParameter(size_t index) const { return this->operator[](index); }
    size_t getParameterCount() const { return getSize(); }

    using array::begin;
    using array::end;

  protected:
    virtual void onSampleRateChanged() {}
    analog_t dt() const { return deltaTime; }
  };

  template<typename T>
  class unitGenerator : public unit, public generator<T>
  {
  protected:
    template<size_t N>
    explicit unitGenerator(init<N>& init, analog_t sampleRate = 1) : unit(init, sampleRate), generator<T>() {}
    unitGenerator(const char_t* name, parameter* params, size_t paramsCount, analog_t sampleRate = 1) : unit(name, params, paramsCount, sampleRate), generator<T>() {}
  };

  template<typename T>
  class unitProcessor : public unit, public processor<T>
  {
  protected:
    template<size_t N>
    explicit unitProcessor(init<N>& init, analog_t sampleRate = 1) : unit(init, sampleRate), processor<T>() {}
    unitProcessor(const char_t* name, parameter* params, size_t paramsCount, analog_t sampleRate = 1) : unit(name, params, paramsCount, sampleRate), processor<T>() {}
  };

  namespace interpolation
  {
    template<typename T>
    struct nearest
    {
      T operator()(const T* buffer, analog_t fracIdx);
    };

    template<typename T>
    struct linear
    {
      T operator()(const T* buffer, analog_t fracIdx);
    };

    template<typename T>
    struct cubic
    {
      T operator()(const T* buffer, analog_t fracIdx);
    };
    
    template<typename T, typename I = linear<T>>
    T sample(const T* buffer, analog_t fracIdx)
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
    T abs(T val) { return ::abs(val); }
    
    template<typename T>
    T constrain(T val, T low, T high) { return val < low ? low : val > high ? high : val; }
    
    template<typename T>
    T epsilon() { return std::numeric_limits<T>::epsilon(); }

    template<typename T>
    T exp(T v) { return ::exp(v); }

    template<typename T>
    T exp2(T v) { return ::exp2(v); }

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

    template<typename T, size_t N>
    frame::channels<T, N> round(frame::channels<T, N> x)
    {
      frame::channels<T, N> result;
      for (size_t i = 0; i < N; i++)
      {
        result[i] = round(x[i]);
      }
      return result;
    }

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
    binary_t isNan(T n) { return isnan(n); }
  }

  namespace random
  {
    // implements the xorshifter algorithm
    static constexpr uint32_t U32_MAX = UINT32_MAX;
    static uint32_t ru32Seed = 33641;
    inline void su32(uint32_t seed) {ru32Seed = seed; }
    inline uint32_t u32()
    {
      ru32Seed ^= ru32Seed << 13; ru32Seed ^= ru32Seed >> 17; ru32Seed ^= ru32Seed << 5;
      return ru32Seed;
    }

    template<typename T>
    T range(T low, T high)
    {
      static constexpr analog_t SCALE = 1/4294967296.0;
      analog_t r = static_cast<analog_t>(u32()) * SCALE;
      return low + r*(high-low);
    }
  }
  
  class gain
  {
    analog_t value;
    explicit gain(analog_t v) : value(v) {}
  public:
    static analog_t decibelsToScale(analog_t db) { return math::exp10(db*static_cast<analog_t>(0.05));}
    static analog_t scaleToDecibels(analog_t scale) { return math::log10(scale)*static_cast<analog_t>(20.0);}
    static gain fromScale(analog_t scale) { return gain(scale); }
    static gain fromDecibels(analog_t dB) { return gain(decibelsToScale(dB)); }

    // implement casting operators so that gain can be used as a parameter type.
    explicit operator binary_t() const { return value > 0; }
    explicit operator digital_t() const { return static_cast<digital_t>(toDecibels()); }
    explicit operator analog_t() const { return value; }

    analog_t toScale() const { return value; }
    analog_t toDecibels() const { return scaleToDecibels(value); }
  };
  
  namespace easing
  {
    struct linear
    {
      analog_t operator()(analog_t t) const { return t; }
    };

    namespace quad
    {
      struct in { analog_t operator()(analog_t t) const { return t * t; } };
      struct out { analog_t operator()(analog_t t) const {  return 1.0 - (1.0 - t) * (1.0 - t); } };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
      struct inOut { analog_t operator()(analog_t t) const { return t < 0.5 ? 2*t*t : 1.0 - math::pow<analog_t>(-2*t+2, 2) * 0.5;  } };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
      struct outIn { static out qo; analog_t operator()(analog_t t) const { return t < 0.5 ? qo(2*t) * 0.5f : 1.0 - qo(2*t) * 0.5; } };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
    }

    namespace expo
    {
      struct in { analog_t operator()(analog_t t) const { return t <= math::epsilon<analog_t>() ? 0 : math::pow<analog_t>(2, 10*t-10); } };
      struct out { analog_t operator()(analog_t t) const { return t >= 1.0 - math::epsilon<analog_t>() ? 1.0 : 1.0 - math::pow<analog_t>(2, -10*t);} };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
      struct inOut
      {
        analog_t operator()(analog_t t) const
        {
          return t <= math::epsilon<analog_t>() ? 0.0
          : t >= 1.0 - math::epsilon<analog_t>() ? 1.0  // NOLINT(cppcoreguidelines-narrowing-conversions)
          : t < 0.5 ? math::pow<analog_t>(2.0, 20*t-10)*0.5  // NOLINT(clang-diagnostic-implicit-float-conversion)
          : 2.0 - math::pow<analog_t>(2.0, -20*t+10)*0.5;  // NOLINT(clang-diagnostic-implicit-float-conversion)
        }
      };
    }

    // easing is first parameter so T can be deduced
    template<typename E, typename T>
    T interp(T begin, T end, analog_t t)
    {
      static E ease;
      return (end-begin) * ease(t) + begin;
    }

    template<typename T>
    T lerp(T begin, T end, analog_t t) { return interp<linear, T>(begin, end, t); }
  }
  
  namespace mixing
  {
    template<typename T, typename E = easing::linear>
    T crossfade(T a, T b, analog_t f)
    {
      static E ease;
      return a*ease(1.0 - f) + b*ease(f);
    }
    
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
          static constexpr size_t coeffSz = 5;
          array<T> src(coeff, coeffSz);
          for (size_t i = 1; i < stages; i++)
          {
            array<T> dst(coeff + coeffSz*i, coeffSz);
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
    static constexpr analog_t B2F = 1.0 / 60.0;  // NOLINT(clang-diagnostic-implicit-float-conversion)
    static constexpr analog_t F2B = 60;

    // can be used by units an indication for how to treat changes in time values (see: delay & freeze)
    enum class mode : uint8_t
    {
      snap, // use duration value directly
      slew, // smooth duration when it changes
      fade, // crossfade to new duration across a block of samples
    };

    analog_t samples; // analog so we can express subsample periods.

    // conversions for parameter
    explicit duration(binary_t b) : samples(b) {}
    explicit duration(size_t s): samples(static_cast<analog_t>(s)) {}
    explicit duration(digital_t i) : samples(static_cast<analog_t>(i)) {}
    explicit duration(analog_t a) : samples(a) {}
    explicit operator binary_t() const { return math::abs(samples) >= math::epsilon<analog_t>(); }
    explicit operator digital_t() const { return static_cast<digital_t>(samples); }
    explicit operator analog_t() const { return samples;}

    static duration fromBpm(analog_t bpm, analog_t sampleRate) { return duration(sampleRate/(bpm*B2F)); }
    static duration fromSeconds(analog_t seconds, analog_t sampleRate) { return duration(sampleRate*seconds); }
    analog_t toBpm(analog_t sampleRate) const { return F2B*(sampleRate/samples); }
    analog_t toSeconds(analog_t sampleRate) const { return samples/sampleRate; }
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
    analog_t sr;

    clockable(analog_t sampleRate, period_t samplePeriodMin, period_t samplePeriodMax, analog_t bpm = 60)
    : tempo(duration::fromBpm(bpm, sampleRate)), periodMin(samplePeriodMin), periodMax(samplePeriodMax)
    , ticks(0), sr(sampleRate)
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
    clockable& operator=(clockable&&) = default;
    // users should call clock at the beginning of every clock pulse
    void clock();
    void clock(period_t sampleDelay);

    analog_t getBpm() const { return tempo.toBpm(sr); }
    // length of one clock pulse in samples
    analog_t getPeriod() const { return tempo.samples;}
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

    virtual T evaluate(analog_t phase) const = 0;  // NOLINT(portability-template-virtual-member-function)
  };

  namespace waves
  {
    template<typename T = double>
    class sine final : public waveform<T>
    {
    public:
      sine() = default;
      T evaluate(analog_t phase) const override { return math::sin(math::twoPi<T>() * phase); }
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
    T evaluate(analog_t phase) const override;
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
    analog_t step;
    T nz[2];

    using unit::dt;

  public:
    explicit noiseGenerator(analog_t sampleRate, noiseTint withTint = noiseTint::white)
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
    analog_t t;

    using unit::dt;

  public:
    explicit ramp(analog_t sampleRate, analog_t durationInSeconds = 0, T fromValue = T(0), T toValue = T(0))
    : unitGenerator<T>(init, sampleRate), mFrom(fromValue), mTo(toValue), t(0)
    {
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

    binary_t isActive() const { return !eor().template read<binary_t>(); }
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
      analog_t time;
      parameter& aw() { return init.params[2]; }
      parameter& ew() { return init.params[3]; }
      
      template<typename E>
      T step();

    public:
      explicit stage(analog_t sampleRate) : unitGenerator<T>(init, sampleRate), mTarget(0), begin(0) { reset(); }

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
    envelope(stage* stages, size_t stageCount, analog_t sampleRate) : unitGenerator<T>(init, sampleRate)
    , stages(stages, stageCount), stageIdx(0), final(sampleRate) {}
    
  public:
    stage& getStage(size_t idx) { return idx == stages.getSize() ? final : stages[idx]; }
    const stage& getStage(size_t idx) const { return idx == stages.getSize() ? final : stages[idx]; }
    size_t getStageCount() const { return stages.size() + 1; }

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
    void startStage(size_t idx, T fromValue) { getStage(idx).start(fromValue); stageIdx = idx; }
    
    // by default, stages advance automatically when their eos goes high.
    // subclasses can override this behavior per stage
    // to enable advancing to the next stage before it is finished,
    // or holding a stage for some period of time.
    virtual binary_t shouldAdvance(size_t currentStageIdx) { return getStage(currentStageIdx).eos().template read<bool>(); }
  };

  template<typename T>
  class ad : public envelope<T>
  {
    typename envelope<T>::stage attackStage;

  public:
    ad(analog_t attackDuration, analog_t decayDuration, analog_t sampleRate) : envelope<T>(&attackStage, 1, sampleRate), attackStage(sampleRate)
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
    binary_t gateOn;
    
  public:
    asr(analog_t attackDuration, analog_t decayDuration, analog_t sampleRate, T triggerThreshold = T(0))
    : ad<T>(attackDuration, decayDuration, sampleRate), triggerThreshold(triggerThreshold), gateOn(false) {}

    using ad<T>::attack;
    typename envelope<T>::stage& release() { return envelope<T>::finalStage(); }
    using ad<T>::eoc; 

    void gate(T value);
    void gate(binary_t on) { gate(on ? T(1) : T(0));}
    void trigger() override { attack().target() << T(1); ad<T>::trigger(); gateOn = false; }
    using envelope<T>::generate;

  protected:
    binary_t shouldAdvance(size_t currentStageIdx) override { return ad<T>::shouldAdvance(currentStageIdx) && !gateOn; }
  };

  template<typename T>
  class adsr : public envelope<T>
  {
    typename envelope<T>::stage attackStage;
    typename envelope<T>::stage decayStage;
    binary_t gateOn;

  public:
    adsr(analog_t attackDuration, analog_t decayDuration, analog_t sustainLevel, analog_t releaseDuration, analog_t sampleRate)
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
    void gate(binary_t on)
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
    binary_t shouldAdvance(size_t currentStageIdx) override
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
    slew(analog_t sampleRate, analog_t riseRate, analog_t fallRate, T initialValue = T(0), T epsilon = math::epsilon<T>()*1000)
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
    explicit smoother(analog_t smoothingDegree = 0.9, T initialValue = T(0)) : unitProcessor<T>(init), mValue(initialValue)  // NOLINT(clang-diagnostic-implicit-float-conversion)
    { degree() << smoothingDegree; }

    parameter& degree() { return init.params[0]; }
    const parameter& value() const { return init.params[1]; }

    T process(const T& v) override { return mValue = easing::lerp(v, mValue, math::constrain<analog_t>(*degree(), 0.0, 1.0)); }
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
    analog_t phase;

    using unit::dt;

  public:
    explicit oscil(analog_t sampleRate, waveform<T>& wave, analog_t freqInHz = 440)
    : unitGenerator<T>(init, sampleRate), wave(&wave), phase(0)
    { fHz() << freqInHz; }

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
    T read(size_t sampleDelay) const;

    // reads behind the write head with a fractional sampleDelay and given interpolation
    template<typename I>
    T read(analog_t sampleDelay) const;

    // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
    T evaluate(analog_t phase) const override;
  };
  
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
    analog_t mDelayInSamples;

    using unit::dt;
    using unit::getSampleRate;

  public:
    delay(array<T> buffer, analog_t sampleRate, analog_t delayInSeconds = 0, analog_t feedbackAmount = 0)
    : unitProcessor<T>(init, sampleRate), mTime(duration::fromSeconds(delayInSeconds, sampleRate))
    , buffer(buffer.getData(), buffer.getSize())
    {
      mDelayInSamples = mTime.samples;
      feedback() << feedbackAmount;
    }

    delayline<T>& getbuffer() { return buffer; }
    const delayline<T>& getBuffer() const { return buffer; }

    /// delay time expressed as vessl::time (i.e. samples), can be set using an analog_t
    parameter& time() { return init.params[0]; }
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    parameter& feedback() { return init.params[1]; }
    
    T process(const T& in) override;

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<duration::mode TimeMode = duration::mode::slew>
    void process(array<T> input, array<T> output);
  };

  template<typename T>
  class follow : public unitProcessor<T>
  {
    unit::init<1> init = {
      "envelope follower",
      {
        parameter("response time", parameter::type::analog),
      }
    };

  public:
    follow(array<T> window, analog_t sampleRate, analog_t responseTimeInSeconds)
      : unitProcessor<T>(init, sampleRate), writer(window), window(window)
      , delta(math::exp(-1.0 / (sampleRate*responseTimeInSeconds))), previous(0), current(0)
    {
      
    }

    T process(const T& in) override
    {
      writer.write(in);
      if (writer.isFull())
      {
        previous = current;
        current = T(0);
        typename
        array<T>::reader r = window.getReader();
        while (r)
        {
          current *= delta;
          current += (1.0 - delta)*math::abs(r.read());
        }
        writer = window.getWriter();
      }

      analog_t t = 1.0 - static_cast<analog_t>(writer.available()) / static_cast<analog_t>(window.getSize());
      return previous + (current - previous) * t;
    }

    using processor<T>::process;

    typename
    array<T>::writer writer;
    array<T> window;
    analog_t delta;
    T previous;
    T current;
  };
  
  // when used as a processor, will write incoming audio to the delayline
  // and output the incoming signal if freeze is not engaged.
  // when freeze is engaged, it will ignore the incoming signal
  // and generate audio from previously recorded input.
  //
  // when used as a generator, it will always generate using the
  // contents of the delayline, enabling it to be used to "freeze"
  // audio that has been recorded elsewhere (e.g. by a delay).
  template<typename T, typename I = interpolation::linear<T>>
  class freeze : public unit, public processor<T>, public generator<T>
  {
    delayline<T> buffer;
    duration mSize;
    analog_t mPhase;
    smoother<analog_t> crossfade; // used to crossfade between the incoming signal and the freeze signal when enabled changes.
    analog_t freezeDelay, freezeSize, readRate;

    init<4> init = {
      "freeze",
      {
        parameter("enabled", parameter::type::binary),
        // end of the freeze loop in samples relative to the most recently recorded sample
        parameter("position", parameter::type::analog),
        // size of the freeze loop as a duration (samples).
        // the beginning of the freeze loop, when played forward,
        // will be position + size.
        parameter("size", &mSize),
        // rate of playback when enabled, can be negative to play in reverse
        parameter("rate", parameter::type::analog)
      }
    };

  public:
    freeze(array<T> buffer, analog_t sampleRate) : unit(init, sampleRate), buffer(buffer.getData(), buffer.getSize())
    , mSize(buffer.getSize()-1), mPhase(0), crossfade(0.75f), freezeDelay(0), freezeSize(mSize.samples), readRate(1)
    { rate() << 1.0; }

    delayline<T>& getBuffer() { return buffer; }
    const delayline<T>& getBuffer() const { return buffer; }
    
    parameter& enabled() { return init.params[0]; }
    parameter& position() { return init.params[1]; }
    parameter& size() { return init.params[2]; }
    parameter& rate() { return init.params[3]; }
    // should this be a parameter? there's not much gained by it.
    analog_t phase() const { return mPhase; }
    // reset the phase to zero (argument for making it a parameter?)
    void reset() { mPhase = 0; }

    T generate() override
    {
      freezeDelay = *position();
      freezeSize  = mSize.samples;
      analog_t sampleDelay = freezeDelay + (1.0-mPhase)*freezeSize;
      mPhase = math::wrap01(mPhase + *rate() / freezeSize);
      return buffer.template read<I>(sampleDelay);
    }
    
    T process(const T& in) override
    {
      binary_t isEnabled = static_cast<binary_t>(enabled());
      analog_t wetLevel = crossfade.process(isEnabled ? 1.0 : 0.0);
      T wet = generate();
      if (!isEnabled)
      {
        buffer.write(in);
      }
      return wetLevel*wet + (1.0 - wetLevel)*in;
    }

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<duration::mode TimeMode = duration::mode::slew>
    void generate(array<T> output) { procgen<TimeMode, false>(output, output); }

    template<duration::mode TimeMode = duration::mode::slew>
    void process(array<T> input, array<T> output) { procgen<TimeMode, true>(input, output); }

  private:
    // shared routine used by templated processing and generation methods
    template<duration::mode TimeMode, bool UseInput>
    void procgen(array<T> input, array<T> output)
    {
      typename array<T>::reader r(input);
      typename array<T>::writer w(output);

      if (TimeMode == duration::mode::snap)
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

      if (TimeMode == duration::mode::slew)
      {
        analog_t fd = *position();
        analog_t fs = mSize.samples;
        analog_t rt = *rate();
        analog_t st = dt();
        while (w)
        {
          freezeDelay = easing::lerp(freezeDelay, fd, st*20);
          freezeSize = easing::lerp(freezeSize, fs, st*20);
          readRate = easing::lerp(readRate, rt, st*10);
          analog_t sampleDelay = freezeDelay + (1.0 - mPhase)*freezeSize;
          mPhase = math::wrap01(mPhase + readRate/freezeSize);
          T wet = buffer.template read<I>(sampleDelay);

          if (UseInput)
          {
            binary_t isEnabled = enabled().readBinary();
            analog_t wetLevel = crossfade.process(isEnabled ? 1.0 : 0.0);
            T in = r.read();
            if (!isEnabled)
            {
              buffer.write(in);
            }
            w << wet*wetLevel + in*(1.0 - wetLevel);
          }
          else
          {
            w << wet;
          }
        }
      }

      if (TimeMode == duration::mode::fade)
      {
        analog_t fd0 = freezeDelay,   fd1 = *position();
        analog_t fs0 = freezeSize,    fs1 = mSize.samples;
        analog_t fade = 0, fadeInc = 1.0 / output.getSize();
        analog_t r0 = readRate, r1 = *rate();
        analog_t p0 = mPhase, dp0 = r0/fs0, dp1 = r1/fs1;
        while (w)
        {
          analog_t sd0 = fd0 + fs0*(1.0-p0);
          analog_t sd1 = fd1 + fs1*(1.0-mPhase);
          T wet = (1.0 - fade)*buffer.template read<I>(sd0) + fade*buffer.template read<I>(sd1);
          mPhase = math::wrap01(mPhase + dp1);
          p0 = math::wrap01(p0 + dp0);
          fade += fadeInc;

          if (UseInput)
          {
            binary_t isEnabled = enabled().readBinary();
            analog_t wetLevel = crossfade.process(isEnabled ? 1.0 : 0.0);
            T in = r.read();
            if (!isEnabled)
            {
              buffer.write(in);
            }
            w << wetLevel*wet + (1.0 - wetLevel)*in;
          }
          else
          {
            w << wet;
          }
        }
        freezeDelay = fd1;
        freezeSize = fs1;
        readRate = r1;
      }
    }
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
    analog_t piosr;
    
    using unit::dt;

  public:
    filter(analog_t sampleRate, analog_t cutoffInHz, analog_t kyu = 0) : unitProcessor<T>(init, sampleRate)
    , piosr(math::pi<analog_t>()/sampleRate)
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
    void onSampleRateChanged() override { piosr = math::pi<analog_t>() * dt(); }
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
    analog_t rateAlpha;

    using unit::dt;

  public:
    bitcrush(analog_t sampleRate, analog_t bitRate, analog_t bitDepth = MaxBits)
      : unitProcessor<T>(init, sampleRate), prevInput(0), currSample(0), rateAlpha(0)
    { rate() << bitRate; depth() << bitDepth; }

    parameter& rate() { return init.params[0]; }
    parameter& depth() { return init.params[1]; }
    parameter& mangle() { return init.params[2]; }

    T process(const T& in) override;

    using unitProcessor<T>::process;
  };

  // A simple peak limiter adapted from pinchenettes/stmlib via DaisySP
  // @todo how best to handle this for non-floating point types (e.g. frame::stereo::analog_t)
  template<typename T>
  class limiter : public unitProcessor<T>
  {
    gain mPreGain;
    unit::init<1> init = {
      "limiter",
      {
        parameter("pre-gain", &mPreGain),
      }
    };
    T mPeak;

  public:
    explicit limiter(gain preGain = gain::fromDecibels(0))
      : unitProcessor<T>(init, 1), mPreGain(preGain), mPeak(0.5)
    { }

    parameter& preGain() { return init.params[0]; }

    T process(const T& in) override;

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
  array<T> array<T>::add(array other, array dest) const
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
  array<T> array<T>::subtract(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "arrays are have different lengths or destination is too small");
    reader a(data, size);
    reader b(other.data, other.size);
    writer c(dest.data, dest.size);
    while (a)
    {
      c << a.read() - b.read();
    }
    return dest;
  }

  template<typename T>
  array<T> array<T>::scale(analog_t value, array dest) const
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
  array<T> array<T>::multiply(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "arrays are have different lengths or destination is too small");
    reader a(data, size);
    reader b(other.data, other.size);
    writer c(dest.data, dest.size);
    while (a)
    {
      c << a.read() * b.read();
    }
    return dest;
  }

  template<typename T>
  frame::channels<T, 1> frame::channels<T, 1>::toMono() const { return channels(samples[0]); }

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
        binary_t state = pv.b & 1;
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
      case type::binary: pv.b = static_cast<binary_t>(value); break;
      case type::digital: pv.i = static_cast<digital_t>(value); break;
      case type::analog: pv.a = static_cast<analog_t>(value); break;
      case type::user: if (pv.u) { *static_cast<T*>(pv.u) = value; } break;
    }
    return *this;
  }
  
  inline void unit::setSampleRate(analog_t sr)
  {
    sampleRate = sr;
    deltaTime = 1.0f / sr;
    onSampleRateChanged();
  }

  inline void clockable::clock()
  {
    tempo.samples = static_cast<analog_t>(math::constrain(ticks, periodMin, periodMax));
    ticks = 0;
    tock(0);
  }

  inline void clockable::clock(period_t sampleDelay)
  {
    tempo.samples = static_cast<analog_t>(math::constrain(ticks + sampleDelay, periodMin, periodMax));
    ticks = 0;
    tock(sampleDelay);
  }

  template<typename T>
  T interpolation::nearest<T>::operator()(const T* buffer, analog_t fracIdx)
  {
    return buffer[static_cast<size_t>(math::round(fracIdx))];
  }

  template<typename T>
  T interpolation::linear<T>::operator()(const T* buffer, analog_t fracIdx)
  {
    analog_t idx;
    analog_t frac = math::mod(fracIdx, &idx);
    size_t x0 = static_cast<size_t>(idx);
    return buffer[x0] + (buffer[x0 + 1] - buffer[x0]) * frac;
  }

  template<typename T>
  T interpolation::cubic<T>::operator()(const T* buffer, analog_t fracIdx)
  {
    static constexpr analog_t DIV6 = (1. / 6.);
    static constexpr analog_t DIV2 = (0.5);

    analog_t idx;
    analog_t f = math::mod(fracIdx, &idx);
    analog_t fm1 = f - 1.f;
    analog_t fm2 = f - 2.f;
    analog_t fp1 = f + 1.f;
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
    analog_t phase = 0;
    analog_t step = 1.0 / N;
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
  T wavetable<T, N, I>::evaluate(analog_t phase) const
  {
    analog_t idx = phase * N;
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
        return 2 * random::range<T>(0, 1) - 1;
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
        analog_t alpha = dt / (dt + RC);
        T white = 2 * random::range(0, 1) - 1;
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
        static uint32_t whiteValues[6] = {
          random::u32() % (RANGE / 6), random::u32() % (RANGE / 6), random::u32() % (RANGE / 6),
          random::u32() % (RANGE / 6), random::u32() % (RANGE / 6), random::u32() % (RANGE / 6)
        };

        int lastKey = key;
        T sum = 0;
        key = key == MAX_KEY ? 0 : ++key;

        int diff = lastKey ^ key;
        for (int i = 0; i < 6; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            whiteValues[i] = random::u32() % (RANGE / 6);
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
    analog_t alpha = step / dt();
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
    if (*duration() > math::epsilon<analog_t>())
    {
      t = 0;
      eorw() << false;
    }
    else
    {
      t = 1;
      eorw() << true;
    }
  }

  template<typename T>
  T ramp<T>::generate()
  {
    analog_t lt = t;
    if (isActive())
    {
      analog_t dinv = 1.0 / *duration();
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
    analog_t dt = unit::dt();
    analog_t s = unit::dt() / math::max<analog_t>(*duration(), dt);
    time += s;
    analog_t t = math::constrain<analog_t>(time, 0.0, 1.0);
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
    binary_t valueOn = value > triggerThreshold;
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
    binary_t isRise = v > value+epsilon;
    binary_t isFall = v < value-epsilon;
    analog_t rate = isRise ? *rise() : isFall ? -1.0f*fall() : (v-value)*getSampleRate();
    value += rate*dt();
    rw() << isRise;
    fw() << isFall;
    return value;
  }

  template<typename T>
  T oscil<T>::generate()
  {
    T val = wave->evaluate(phase + *pm());
    phase += (*fHz() * math::exp2(*fmExp()) + *fmLin()) * dt();
    phase = math::wrap01(phase);
    return val;
  }

  template<typename T>
  T delayline<T>::read(size_t sampleDelay) const
  {
    assert(sampleDelay < getSize());
    sampleDelay = getSize() - 1 - sampleDelay;
    size_t idx = getWriteIndex() + sampleDelay;
    return getData()[idx % getSize()];
  }

  template<typename T>
  template<typename I>
  T delayline<T>::read(analog_t sampleDelay) const
  {
    assert(sampleDelay >= 0 && sampleDelay < getSize());
    analog_t fSize = static_cast<analog_t>(getSize());
    sampleDelay = fSize - 1 - sampleDelay;
    analog_t fidx = static_cast<analog_t>(getWriteIndex()) + sampleDelay;
    analog_t idx;
    analog_t f = math::mod(fidx, &idx);
    size_t x0 = static_cast<size_t>(idx) % getSize();
    size_t x1 = (x0 + 1) % getSize();
    size_t x2 = (x0 + 2) % getSize();
    const T* data = getData();
    T s[3] = { data[x0], data[x1], data[x2] };
    return interpolation::sample<T, I>(s, f);
  }

  template<typename T>
  T delayline<T>::evaluate(analog_t phase) const
  {
    analog_t fSize = static_cast<analog_t>(getSize());
    phase = math::wrap<analog_t>(phase, -1.0, 1.0);
    analog_t sampleDelay = phase > 0 ? (1.0 - phase) * fSize : -phase * fSize;
    return read<interpolation::linear<T>>(sampleDelay);
  }

  template<typename T, typename I>
  T delay<T, I>::process(const T& in)
  {
    mDelayInSamples = mTime.samples;
    // delay time in samples
    analog_t dts = math::constrain<analog_t>(mDelayInSamples, 0.f, static_cast<analog_t>(buffer.getSize()-1));
    analog_t s = buffer.template read<I>(dts);
    analog_t fbk = math::constrain<analog_t>(*feedback(), -1.0, 1.0);
    buffer.write(in + s * fbk);
    return s;
  }

  template<typename T, typename I>
  template<duration::mode TimeMode>
  void delay<T, I>::process(array<T> input, array<T> output)
  {
    if (TimeMode == duration::mode::snap)
    {
      processor<T>::process(input, output);
    }

    if (TimeMode == duration::mode::slew)
    {
      typename array<T>::reader r = input.getReader();
      typename array<T>::writer w = output.getWriter();
      analog_t dst = dt() * 10.0f;
      while (r && w)
      {
        T in = r.read();
        mDelayInSamples = easing::lerp(mDelayInSamples, mTime.samples, dst);
        // delay time in samples
        analog_t dts = math::constrain<analog_t>(mDelayInSamples, 0.f, static_cast<analog_t>(buffer.getSize()-1));
        analog_t wet = buffer.template read<I>(dts);
        analog_t fbk = math::constrain<analog_t>(*feedback(), -1.0, 1.0);
        buffer.write(in + wet*fbk);
        w << wet;
      }
    }
    
    if (TimeMode == duration::mode::fade)
    {
      typename array<T>::reader r = input.getReader();
      typename array<T>::writer w = output.getWriter();
      
      analog_t fade = 0;
      analog_t fadeInc = 1.0f / input.getSize();
      // smooth time parameter to prevent crunchiness when it is noisy or changes by large amounts
      analog_t targetSampleDelay = mTime.samples;
      // delay time in samples
      analog_t fts = math::constrain<analog_t>(mDelayInSamples, 0.0, static_cast<analog_t>(buffer.getSize()-1));
      analog_t tts = math::constrain<analog_t>(targetSampleDelay, 0.0, static_cast<analog_t>(buffer.getSize()-1));
      analog_t fbk = math::constrain<analog_t>(*feedback(), -1.0, 1.0);
      
      while (r && w)
      {
        T in = r.read();
        T wet = (1.0f - fade) * buffer.template read<I>(fts) + fade * buffer.template read<I>(tts);
        buffer.write(in + wet*fbk);
        w << wet;
        fade += fadeInc;
      }
        
      mDelayInSamples = targetSampleDelay;
    }
  }

  template<typename T, uint32_t MaxBits>
  T bitcrush<T, MaxBits>::process(const T& in)
  {
    rateAlpha += math::max(1.0f, *rate())*dt();
    if (rateAlpha >= 1)
    {
      rateAlpha -= 1;
      currSample = easing::lerp(prevInput, in, rateAlpha);
    }

    analog_t bd = math::constrain<analog_t>(*depth(), 2.0, MaxBits);
    analog_t scalar = math::pow<analog_t>(2.0, bd) - 1;
    T val = math::round(currSample*scalar);
    if (mangle().template read<bool>())
    {
      val = val ^ math::round(prevInput*scalar);
    }
    prevInput = in;
    return val * (1.0 / scalar);
  }

  template<typename T>
  T limiter<T>::process(const T& in)
  {
    T pre = in*mPreGain.toScale();
    T peak = math::abs(pre);
    T error = peak - mPeak;
    mPeak += (error > T(0)) ? T(0.05) : T(0.00002) * error;
    analog_t gain = (mPeak <= 1.0 ? 1.0 : 1.0 / mPeak);
    // DaisySP returns this, which sounds better for how I typically use this.
    return saturation::softlimit(pre*gain*0.7);
    // stmlib returns this, which clips more easily, but is faster.
    //return pre*gain*0.8;
  }

  template<>
  inline float math::wrap01(float v) { float i; float f = modf(v, &i); return f < 0 ? 1.0f - f : f; }

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
  array<float> array<float>::add(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_add_f32(data, other.data, dest.data, size);
    return dest;
  }

  template<>
  array<float> array<float>::subtract(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_sub_f32(data, other.data, dest.data, size);
    return dest;
  }

  template<>
  array<float> array<float>::scale(float value, array dest) const
  {
    VASSERT(size <= dest.size, "Destination array is not large enough");
    arm_scale_f32(data, value, dest.data, size);
    return dest;
  }

  template<>
  array<float> array<float>::multiply(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_mult_f32(data, other.data, dest.data, size);
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
