////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2025 Damien Quartz
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
// ReSharper disable once CppUnusedIncludeDirective
#include <cstring>
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

#ifdef exp10
#undef exp10
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
  using size_t = uint64_t;
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
    sink& operator<<(const T& value) { write(value); return *this; }
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
  class list
  {
  public:
    virtual ~list() = default;
    virtual size_t getSize() const = 0;
    binary_t isEmpty() const { return getSize() == 0; }
    T operator[](size_t index) const
    {
      VASSERT(index < getSize(), "Attempt to access a list element with out-of-bounds index."); 
      return elementAt(index);
    }
    
    class iterator
    {
      const list* src = nullptr;
      size_t index = 0;
      iterator(const list& src, size_t idx): src(&src), index(idx) {}
    public:
      static iterator begin(const list& src) { return iterator(src, 0); }
      static iterator end(const list& src) { return iterator(src, src.getSize()); }
      T operator*() const { return (*src)[index]; }
      iterator& operator++() { index++; return *this; } 
      bool operator==(const iterator& it) const { return src == it.src && index == it.index; }
      bool operator!=(const iterator& it) const { return !(*this == it); }
    };
    
  protected:
    virtual T elementAt(size_t index) const = 0;
  };
  
  template<typename T>
  typename list<T>::iterator begin(const list<T>& lst) { return list<T>::iterator::begin(lst); }
  
  template<typename T>
  typename list<T>::iterator end(const list<T>& lst) { return list<T>::iterator::end(lst); }

  template<typename T>
  class array // not a list<T> to keep size to 16 bytes
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

    // ranged-based for support
    T* begin() { return data; }
    T* end() { return data + size; }
    
    T& operator[](size_t index) { return data[index]; }
    const T& operator[](size_t index) const { return data[index]; }

    class reader final : public source<T>
    {
    protected:
      const T* begin;
      const T* head;
      const T* end;

    public:
      explicit reader() : source<T>(), begin(nullptr), head(nullptr), end(nullptr) {}
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
      // block copy the entire contents of reader into this writer.
      // writer must have enough space for the contents of reader.
      // ReSharper disable once CppEnforceOverridingFunctionStyle
      void write(const reader& r);
      size_t available() const { return end - head; }
    };

    reader getReader() const { return reader(*this); }
    writer getWriter() { return writer(*this); }

    // block copy this array to dest, which must be large enough to hold this array.
    void copyTo(array dest);

    void fill(T value);
    
    // adds value to every element in the array, returns dest
    array offset(T value, array dest) const;
    array offset(T value) { return offset(value, *this); }
    
    // element-wise addition of this and other, returns dest
    array add(array other, array dest) const;
    array add(array other) { return add(other, *this); }
    
    // element-wise subtraction of this and other, returns dest
    array subtract(array other, array dest) const;
    array subtract(array other) { return subtract(other, *this); }
    
    // scales every element in this array by value, returns dest
    array scale(analog_t value, array dest) const;
    array scale(analog_t value) { return scale(value, *this); }
    
    // element-wise multiplication of this and other, returns dest
    array multiply(array other, array dest) const;
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
      channels(channels const& other) : array<T>(samples, N) { other.copyTo(*this); }
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
      ~channels() = default;
      
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
      ~channels() = default;
      
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
  
  // for the most part these resolve to standard math calls,
  // template specializations are provided for ARM_CORTEX where applicable.
  namespace math
  {
    template<typename T>
    constexpr T e() { return static_cast<T>(2.71828182845904523536); }
    
    template<typename T>
    constexpr T pi() { return static_cast<T>(3.1415926535897932385); }

    template<typename T>
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
    T exp10(T v) { return ::pow(T(10), v); }
    
    template<typename T>
    T log(T v) { return ::log(v); }

    template<typename T>
    T log10(T v) { return ::log10(v); }
    
    template<typename T>
    T max(const T& a, const T& b) { return a > b ? a : b; }

    template<typename T>
    T min(const T& a, const T& b) { return a < b ? a : b; }

    template<typename T>
    T mod(T v, T* i) { return ::modf(v, i); }

    template<typename T>
    T pow(T x, T y) { return ::pow(x, y); }

    template<typename T>
    T round(T x) { return ::round(x); }
    
    template<typename T>
    T floor(T x) { return ::floor(x); }

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
    T cos(T x) { return ::cos(x); }

    template<typename T>
    T sqrt(T x) { return ::sqrt(x); }
    
    template<typename T>
    T sqrt2() { static T v = sqrt(2); return v; }

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
    
    template<>
    inline float wrap01(float v) { float i; float f = modf(v, &i); return f < 0 ? 1.0f - f : f; }

    template<typename T>
    binary_t isNan(T n) { return isnan(n); }
    
    // xore because xor is a keyword
    template<typename T>
    T xore(const T& a, const T& b)
    {
      return a ^ b;
    }
    
    template<>
    inline analog_t xore(const analog_t& a, const analog_t& b)
    {
      return xore(static_cast<digital_t>(a), static_cast<digital_t>(b));
    }
  }
  
  // stored as decibels (0 = unity gain).
  struct gain
  {
    analog_t db;
    
    gain() : db(0) {}
    explicit gain(analog_t dbv) : db(dbv) {}
    static analog_t decibelsToScale(analog_t db) { return math::exp10(db*static_cast<analog_t>(0.05));}
    static analog_t scaleToDecibels(analog_t scale) { return math::log10(scale)*static_cast<analog_t>(20.0);}
    static gain fromScale(analog_t scale) { return gain(scaleToDecibels(scale)); }
    static gain fromDecibels(analog_t dB) { return gain(dB); }

    // implement casting operators so that gain can be used as a parameter type.
    explicit operator binary_t() const { return db >= 0; }
    explicit operator digital_t() const { return static_cast<digital_t>(db); }
    explicit operator analog_t() const { return db; }

    analog_t toScale() const { return decibelsToScale(db); }
    analog_t toDecibels() const { return db; }
  };
  
  struct duration
  {
    // convert bpm to frequency in Hz
    static constexpr analog_t B2F = 1.0 / 60.0;  // NOLINT(clang-diagnostic-implicit-float-conversion)
    static constexpr analog_t F2B = 60;

    // can be used by units as an indication for how to treat changes in time values (see: delay & freeze)
    enum class mode : uint8_t
    {
      snap, // use duration value directly
      slew, // smooth duration when it changes
      fade, // crossfade to new duration across a block of samples
    };

    analog_t samples; // analog so we can express subsample periods.

    duration() : samples(0) {}
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
  struct procSample
  {
    processor<T>* proc;
    T sample;
    procSample(processor<T>& proc, T sample) : proc(&proc), sample(sample) {}
    T& operator>>(T& out) { out = proc->process(sample); return out; }
  };
  
  template<typename T>
  struct procBlock
  {
    processor<T>* proc;
    array<T> block;
    procBlock(processor<T>& proc, array<T> arr) : proc(&proc), block(arr) {}
    // implies in-place processing of block (but what if it didn't? how do?)
    procBlock operator>>(processor<T>& next) { proc->process(block, block); return procBlock(next, block); }
    array<T> operator>>(array<T> out) { proc->process(block, out); return out; }
  };

  template<typename T>
  procBlock<T> operator>>(array<T> in, processor<T>& proc) { return procBlock<T>(proc, in); }
  
  template<typename T>
  struct procStream : source<T>
  {
    processor<T>* proc;
    source<T>* src;
    procStream(processor<T>& proc, source<T>& src) : proc(&proc), src(&src) {}
    binary_t isEmpty() const override { return src->isEmpty(); }
    T read() override { return proc->process(src->read()); }
    T& operator>>(T& rhs) { rhs = read(); return rhs; }
    procStream operator>>(processor<T>& proc) { return procStream(proc, *this); }
    sink<T>& operator>>(sink<T>& out) { proc->process(*src, out); return out; }
    array<T> operator>>(array<T> out) { typename array<T>::writer w = out.getWriter(); proc->process(*src, w); return out; }
  };
  
  template<typename T>
  procBlock<T> operator>>(const T& in, processor<T>& proc) { return procBlock<T>(proc, frame::channels<T,1>(in)); }

  template<typename T>
  procStream<T> operator>>(source<T>& in, processor<T>& proc) { return procStream<T>(proc, in); }

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
    enum class valuetype : uint8_t
    {
      none = 0,
      binary = 1, // on/off (binary_t)
      digital = 2, // integral values (digital_t)
      analog = 3, // floating point values (analog_t)
      // space for more built-ins

      // user provided type, stored as a void*, must be convertible to all other parameter types.
      // even if the conversion is meaningless.
      user = UINT8_MAX 
    };
    
    typedef uint32_t id_t;
    
    struct desc
    {
      const char_t* name;
      id_t          id;
      valuetype     type;
      
      static desc empty() { return {"", 0, valuetype::none}; }
    };
    
    template<size_t N>
    struct desclist
    {
      static constexpr size_t size = N;
      desc descs[N];
      
      constexpr desc operator[](id_t id) const
      {
        for (size_t i = 0; i < size; ++i)
        {
          if (descs[i].id == id)
          {
            return descs[i];
          }
        }
        return desc::empty();
      }
    };
    
    template<typename T>
    struct data
    {
      T value = T(0);
      static constexpr auto type = valuetype::user;
    };
    
    template<typename T>
    parameter(const desc& inDesc, const data<T>& inData) : description(inDesc), pdata(&const_cast<data<T>&>(inData).value) {}
    
    const desc& getDescription() const { return description; }
    
    template<typename T>
    T read() const
    {
      return *static_cast<const T*>(pdata);
    }

    binary_t  readBinary()  const { return read<binary_t>(); }
    digital_t readDigital() const { return read<digital_t>(); }
    analog_t  readAnalog()  const { return read<analog_t>(); }
    
    explicit operator binary_t()  const { return read<binary_t>(); }
    explicit operator digital_t() const { return read<digital_t>(); }
    explicit operator analog_t()  const { return read<analog_t>(); }

    // static-cast T to parameter type before assign
    template<typename T>
    parameter& write(const T& value)
    {
      switch (description.type)
      {
        case valuetype::none: break;
        case valuetype::binary: *static_cast<binary_t*>(pdata)   = static_cast<binary_t>(value); break;
        case valuetype::digital: *static_cast<digital_t*>(pdata) = static_cast<digital_t>(value); break;
        case valuetype::analog: *static_cast<analog_t*>(pdata)   = static_cast<analog_t>(value); break;
        // attempt to cast data to T, might work?
        case valuetype::user: *static_cast<T*>(pdata) = value; break;
      }
      return *this;
    }
    
    template<typename T>
    parameter& operator=(const T& value)
    {
      write(value);
      return *this;
    }
    
    parameter& operator=(const parameter& rhs)
    {
      if (&rhs != this)
      {
        switch (description.type)
        {
        case valuetype::none: break;
        case valuetype::binary: *static_cast<binary_t*>(pdata)   = rhs.readBinary(); break;
        case valuetype::digital: *static_cast<digital_t*>(pdata) = rhs.readDigital(); break;
        case valuetype::analog: *static_cast<analog_t*>(pdata)   = rhs.readAnalog(); break;
        case valuetype::user: 
          switch (rhs.description.type)
          {
        case valuetype::none: break;
        case valuetype::binary: write(rhs.readBinary()); break;
        case valuetype::digital: write(rhs.readDigital()); break;
        case valuetype::analog: write(rhs.readAnalog()); break;
        case valuetype::user: VASSERT(false, "Can't assign user parameter to user parameter with operater="); break;
          }
          break;
        }
      }
      return *this;
    }
    
    static parameter none();
    
  private:
    desc  description;
    void* pdata;
  };
  
  struct parameters : list<parameter>
  {
    // @todo access by ID
  };

  inline parameters::iterator begin(const parameters& lst) { return parameters::iterator::begin(lst); }
  inline parameters::iterator end(const parameters& lst)   { return parameters::iterator::end(lst); }
  
  template<size_t N>
  struct plist : parameters
  {
    static constexpr size_t plsz = N;
    size_t getSize() const override { return N; }
  };
  
  template<>
  struct parameter::data<void*>
  {
    void* value = nullptr;
    static constexpr auto type = valuetype::none;
  };
  
  template<>
  struct parameter::data<analog_t>
  {
    analog_t value;
    static constexpr auto type = valuetype::analog;
  };
  
  template<>
  struct parameter::data<digital_t>
  {
    digital_t value;
    static constexpr auto type = valuetype::digital;
  };
  
  template<>
  struct parameter::data<binary_t>
  {
    binary_t value;
    static constexpr auto type = valuetype::binary;
  };
  
  template<>
  struct parameter::data<gain>
  {
    gain value = gain();
    // @todo add a gain type
    static constexpr auto type = valuetype::user;
  };
  
  template<>
  struct parameter::data<duration>
  {
    duration value = duration();
    // @todo add a duration type
    static constexpr auto type = valuetype::user;
  };
  
  template<typename T>
  struct param : parameter::data<T>
  {
    parameter operator()(const parameter::desc& d) const { return parameter(d, *this); }
  };
  
  typedef param<analog_t>  analog_p;
  typedef param<digital_t> digital_p;
  typedef param<binary_t>  binary_p;
  typedef param<gain>      gain_p;
  typedef param<duration>  duration_p;

  inline parameter parameter::none() { data<void*> v; return parameter(desc::empty(), v); }
  
  template<typename T>
  procSample<T> operator>>(const parameter& p, processor<T>& proc) { return procSample<T>(proc, p.read<T>()); }
  
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
  T operator+(const T& v, const parameter& p) { return v + p.read<T>(); }
  
  template<typename T>
  T operator-(const parameter& p, const T& v) { return p.read<T>() - v; }
  
  template<typename T>
  T operator-(const T& v, const parameter& p) { return v - p.read<T>(); }
  
  template<typename T>
  T operator*(const parameter& p, const T& v) { return p.read<T>() * v; }
  
  template<typename T>
  T operator*(const T& v, const parameter& p) { return v * p.read<T>(); }
  
  template<typename T>
  T operator/(const parameter& p, const T& v) { return p.read<T>() / v; }
  
  template<typename T>
  T operator/(const T& v, const parameter& p) { return v / p.read<T>(); }

  class unit
  {
  public:
    struct description
    {
      const char_t*                name;
      const parameter::desc*       params;
      size_t                       paramCount;
    };

    virtual ~unit() = default;
    
    // common interface for setting sample rate, for those units that might need it.
    virtual void setSampleRate(analog_t sr) {}
    
    // providing a description is optional, but useful.
    virtual description getDescription() const { return { "", nullptr, 0 }; }
    virtual const parameters& getParameters() const = 0;
    
    // @todo get parameter by name / id
  };
  
  inline const parameter::desc* begin(const unit::description& desc) { return desc.params; }
  inline const parameter::desc* end(const unit::description& desc) { return desc.params + desc.paramCount; }

  template<typename T>
  class unitGenerator : public unit, public generator<T>
  {
  protected:
    explicit unitGenerator() : unit(), generator<T>() {}
  };

  template<typename T>
  class unitProcessor : public unit, public processor<T>
  {
  protected:
    explicit unitProcessor() : unit(), processor<T>() {}
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
  
  namespace easing
  {
    struct linear
    {
      analog_t operator()(analog_t t) const { return t; }
    };

    struct smoothstep
    {
      analog_t operator()(analog_t t) const { return t * t * (3.0f - 2.0f * t); }
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
    
    template<typename T>
    T smooth(T value, T target, analog_t degree = 0.9f) { return value*degree + (1.0 - degree)*target; }
    
    template<>
    inline digital_t smooth(digital_t value, digital_t target, analog_t degree) { return (value*degree + target)/(degree+1);  }
  }

  // analog unipolar noise generators that generate values in the range [0,1]
  namespace noise
  {
    struct white
    {
      explicit white(analog_t sampleRate) {}
      analog_t operator()() const { return random::range<analog_t>(0, 1); }
    };

    // Implements the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
    // Would be good to dig into the improvements on the algorithm mentioned later in the article.
    struct pink
    {
    private:
      static constexpr int RANGE = 128;
      static constexpr int COUNT = 6;
      static constexpr int MAX_KEY = 0x1f;

      int key = 0;
      analog_t maxSum = 90;
      uint32_t whiteValues[COUNT] = {
        random::u32() % (RANGE / COUNT), random::u32() % (RANGE / COUNT), random::u32() % (RANGE / COUNT),
        random::u32() % (RANGE / COUNT), random::u32() % (RANGE / COUNT), random::u32() % (RANGE / COUNT)
      };
      
    public:
      explicit pink(analog_t sampleRate) {}
      
      analog_t operator()()
      {
        int lastKey = key;
        analog_t sum = 0;
        key = key == MAX_KEY ? 0 : ++key;

        int diff = lastKey ^ key;
        for (int i = 0; i < COUNT; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            whiteValues[i] = random::u32() % (RANGE / COUNT);
          }
          sum += static_cast<analog_t>(whiteValues[i]);
        }
        maxSum = math::max(sum, maxSum);
        analog_t n = sum / maxSum;
        // vassert(!math::isNan(n) && "pink noise generated nan");
        return n;
      }
    };

    // Brownian noise (i.e. random wander) run thru a DC blocking filter.
    // See: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
    // @todo still a bit crunchy
    struct red
    {
      explicit red(analog_t sampleRate) : r((sampleRate-2.0f)/sampleRate), rc((1.0f - r)*200), x1(0), y1(0) {}
      
      analog_t operator()()
      {
        analog_t white = random::range<analog_t>(-rc, rc);
        analog_t x = y1 + white;
        // only run the filter when we get close to going out of range
        // to compensate for wandering away from 0.
        y1 = x < -0.49 || x > 0.49f ? x - x1 + r*y1 : x;
        x1 = x;
        return y1 + 0.5f;
      }

    private:
      analog_t r;
      analog_t rc;
      analog_t x1, y1;
    };
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
    struct args
    {
      analog_t sr;
      analog_t hz;
      analog_t q;
      gain     g;
      args(analog_t sr, analog_t hz, analog_t q, gain g) : sr(sr), hz(hz), q(q), g(g) {}
      // helper for biquad
      analog_t omega() const { return hz * math::pi<analog_t>() / sr; }
    };
    
    // DC blocking filter, see: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
    template<typename T>
    struct dcblock
    {
      T x1 = T(0), y1 = T(0);
      void process(const T* source, T* dest, size_t blockSize, const args& args)
      {
        // gives us about 0.995 for 44100, which is a pretty good R according to the article above.
        analog_t r = (args.sr - 200.0f) / args.sr;
        while (blockSize--)
        {
          T x = *source++;
          y1 = x - x1 + r*y1;
          x1 = x;
          *dest++ = y1;
        }
      }
    };
    
    // common filter Q values for biquad filters
    namespace q
    {
      template<typename T>
      constexpr T butterworth() { return static_cast<T>(0.70710678118); } // 1/sqrt(2)

      template<typename T>
      constexpr T sallenKey() { return static_cast<T>(0.5); } 

      template<typename T>
      constexpr T bessel() { return static_cast<T>(0.57735026919); } // 1/sqrt(3)
    }
    
    template<typename T>
    struct data
    {
      array<analog_t> coeff;
      array<T> state;

      data(analog_t* coeffData, size_t coeffSize, T* stateData, size_t stateSize)
      : coeff(coeffData, coeffSize), state(stateData, stateSize)
      { coeff.fill(0); state.fill(T(0)); }
      
      size_t getCoeffSize() const { return coeff.getSize(); }
      size_t getStateSize() const { return state.getSize(); }
    };

    // based on https://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
    template<size_t STAGES>
    struct biquad
    {
      static constexpr size_t COEFF_NUM = 5;
      
      template<typename T, size_t STATES>
      struct cascade : data<T>
      {
      private:
        analog_t co[COEFF_NUM*STAGES];
        T st[STATES*STAGES];
      public:
        cascade() : data<T>(co, COEFF_NUM * STAGES, st, STATES * STAGES), co{} {}
      };
      
      template<typename T, class CoGen>
      struct df2T final : cascade<T, 2>
      {
        using cascade<T, 2>::coeff;
        using cascade<T, 2>::state;
        using cascade<T, 2>::getCoeffSize;
        using cascade<T, 2>::getStateSize;
        
        static CoGen cg;
         
        void process(const T* source, T* dest, size_t blockSize, const args& args);
        // ReSharper disable once CppMemberFunctionMayBeStatic
        size_t getStageCount() const { return STAGES; }
      };

      template<typename T>
      static void copy(T* coeff);
      
      // coefficient generators
      struct lpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const; };
      struct hpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const; };
      struct bpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const; };
      struct ntcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const; };
      struct pkcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain g) const; };
      struct lscg { void operator()(analog_t* coeff, analog_t omega, analog_t _, gain g) const; };
      struct hscg { void operator()(analog_t* coeff, analog_t omega, analog_t _, gain g) const; };
      
      // base class for filter types to wrap df2T because we specialize it for ARM.
      template<typename T, class CoGen>
      struct flt
      {
        df2T<T, CoGen> df2;
        void process(const T* source, T* dest, size_t blockSize, const args& args) { df2.process(source, dest, blockSize, args); }
      };

      // filter types for the filter unit generator
      template<typename T>
      struct lowPass final: flt<T, lpcg> {};
      
      template<typename T>
      struct highPass final : flt<T, hpcg> {};
      
      template<typename T>
      struct bandPass final : flt<T, bpcg> {};
        
      template<typename T>
      struct notch final : flt<T, ntcg> {};
        
      template<typename T>
      struct peak final : flt<T, pkcg> {};
        
      template<typename T>
      struct lowShelf final : flt<T, lscg> {};
        
      template<typename T>
      struct highShelf final : flt<T, hscg> {};
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
  
  template<typename T = analog_t>
  struct smoother
  {
    T        value;
    analog_t degree;
    
    explicit smoother(analog_t smoothingDegree = 0.9f, T initialValue = T(0))
      : value(initialValue), degree(smoothingDegree)
    {
    }
    
    explicit operator T() const { return value; }
    
    // so we can use this like OWL's SmoothValue
    T operator=(const T& v)
    {
      analog_t d = math::constrain<analog_t>(degree, 0.0, 1.0);
      return value = easing::smooth(value, v, d);
    }
    
    T operator=(const parameter& p)
    {
      return *this = p.read<T>();
    }
  };
  
  // a waveform that can be evaluated using a normalized phase value
  // implementors should accept negative phase, as well as phase values outside [-1,1]
  template<typename T>
  struct waveform
  {
    using SampleType = T;
    
    waveform() = default;
    virtual ~waveform() = default;
    waveform(const waveform&) = default;
    waveform(waveform&&) = default;
    waveform& operator=(const waveform&) = default;
    waveform& operator=(waveform&&) = default;

    virtual T evaluate(analog_t phase) const = 0;  // NOLINT(portability-template-virtual-member-function)
  };

  namespace waves
  {
    template<typename T = analog_t>
    struct sine final : waveform<T>
    {
      T evaluate(analog_t phase) const override { return math::sin(math::twoPi<T>() * phase); }
    };

    template<typename T = analog_t>
    struct square final : waveform<T>
    {
      analog_t pulseWidth;
      square() : pulseWidth(0.5) {}
      explicit square(analog_t pulseWidth) : pulseWidth(pulseWidth) {}
      T evaluate(analog_t phase) const override { return phase < pulseWidth ? 1 : -1; }
    };

    // same as square, but unipolar
    template<typename T = analog_t>
    struct clock final : waveform<T>
    {
      analog_t pulseWidth;
      clock() : pulseWidth(0.5) {}
      explicit clock(analog_t pulseWidth) : pulseWidth(pulseWidth) {}
      T evaluate(analog_t phase) const override { return phase < pulseWidth ? 1 : 0; }
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
  
  // generates stepped analog noise at rate values per second.
  template<typename T, typename N>
  class noiseGenerator final : public unitGenerator<T>, protected plist<1>
  {
  public:
    explicit noiseGenerator(analog_t sampleRate = 1)
    : unitGenerator<T>(), noise(sampleRate), dt(1.0f/sampleRate), step(0)
    {
      value = noise();
      next = noise();
      params.rate.value = sampleRate;
    }
    
    void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate;}
    const parameters& getParameters() const override { return *this; }

    parameter rate() const { return params.rate({ "rate", 'r', analog_p::type }); }

    // generates stepped noise in the range [0,1] at the given rate
    T generate() override;

    // smooths the stepped noise with the given easing
    template<typename E>
    T generate();

  protected:
    parameter elementAt(size_t index) const override { parameter p[plsz] = { rate() }; return p[index]; }

  private:
    struct
    {
      analog_p rate;
    } params;
    N noise;
    analog_t dt;
    analog_t value;
    analog_t next;
    analog_t step;
  };

  // unit that generates a linear ramp from one value to another over a duration of seconds
  // @todo implement easings above and add that as a template parameter
  template<typename T>
  class ramp final : public unitGenerator<T>, protected plist<4>
  {
  public:
    explicit ramp(analog_t sampleRate, analog_t durationInSeconds = 0, T fromValue = T(0), T toValue = T(0))
    : unitGenerator<T>(), dt(1.0f/sampleRate), t(0)
    {
      params.from.value = fromValue;
      params.to.value = toValue;
      params.duration.value = durationInSeconds;
      params.eor.value = false;
    }
    
    void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate;}
    const parameters& getParameters() const override { return *this; }

    // ins
    parameter from() const { return params.from({ "from", 'f', param<T>::type }); }
    parameter to() const { return params.to({ "to", 't', param<T>::type }); }
    parameter duration() const { return params.duration({ "duration", 'd', analog_p::type }); }

    // outs
    parameter eor() const { return params.eor({ "eor", 'e', binary_p::type }); }
    // could also add t as an out.

    binary_t isActive() const { return !params.eor.value; }
    void trigger();
    T generate() override;
    
  protected:
    parameter elementAt(size_t index) const override { parameter p[plsz] = { from(), to(), duration(), eor() }; return p[index]; }
    
  private:
    struct
    {
      param<T> from, to;
      analog_p duration;
      binary_p eor;
    } params;
    analog_t dt;
    analog_t t;
  };

  // generates an envelope that begins and ends at zero, with some number of stages leading up to a final stage.
  // each stage of an envelope is defined by a target value, a duration to reach it,
  // and whether the envelope should hold the stage value until it is triggered again.
  // @todo figure out how the heck to return a full description and parameter lists for envelopes.
  template<typename T>
  class envelope : public unitGenerator<T>, protected plist<2>
  {
  public:
    class stage final : public unitGenerator<T>, protected plist<5>
    {
    public:
      explicit stage(analog_t sampleRate) : unitGenerator<T>(), begin(0), dt(1.0f/sampleRate) { reset(); }
      
      void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate;}
      const parameters& getParameters() const override { return *this; }

      parameter target() const { return params.target({ "target", 't', analog_p::type }); }
      parameter duration() const { return params.duration({ "duration", 'd', analog_p::type }); }
      parameter active() const { return params.active({ "active", 'a', binary_p::type }); }
      parameter eos() const { return params.eos({ "eos", 'e', binary_p::type }); }
      // current value of the stage
      parameter value() const { return params.output({ "value", 'v', analog_p::type }); }

      void start(T fromValue) { begin = fromValue; params.output.value = fromValue; time = -dt; params.active.value = true; params.eos.value = false; }
      void reset() { params.active.value = false; params.eos.value = false; time = -dt; params.output.value = 0; }
      
      template<typename E>
      T generate() { return active() ? step<E>() : params.output.value; }
      T generate() override { return generate<easing::linear>(); }
      
    protected:
      parameter elementAt(size_t index) const override
      {
        parameter p[plsz] = { target(), duration(), active(), eos(), value() };
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
      } params;
      T begin; // value the stage started with
      // where we are in the stage
      analog_t time;
      analog_t dt;
    };
  protected:
    envelope(stage* stages, size_t stageCount, analog_t sampleRate) : unitGenerator<T>()
    , stages(stages, stageCount), stageIdx(0), final(sampleRate) {}
    
  public:
    void setSampleRate(float sampleRate) override;
    
    stage& getStage(size_t idx) { return idx == stages.getSize() ? final : stages[idx]; }
    const stage& getStage(size_t idx) const { return idx == stages.getSize() ? final : stages[idx]; }
    size_t getStageCount() const { return stages.size() + 1; }

    stage& currentStage() { return getStage(stageIdx); }
    const stage& currentStage() const { return getStage(stageIdx);}
    stage& finalStage() { return final; }
    const stage& finalStage() const { return final; }
    
    const parameters& getParameters() const override { return *this; }

    parameter value() const { return currentStage().value(); }
    parameter eoc() const { return params.eoc({ "eoc", 'e', binary_p::type }); }

    // make this a parameter we check in generate?
    virtual void trigger();

    template<typename E>
    T generate();
    T generate() override { return generate<easing::linear>(); }

  protected:
    parameter elementAt(size_t index) const override { parameter p[plsz] = { value(), eoc() }; return p[index]; }
    void startStage(size_t idx, T fromValue) { getStage(idx).start(fromValue); stageIdx = idx; }
    
    // by default, stages advance automatically when their eos goes high.
    // subclasses can override this behavior per stage
    // to enable advancing to the next stage before it is finished,
    // or holding a stage for some period of time.
    virtual binary_t shouldAdvance(size_t currentStageIdx) { return getStage(currentStageIdx).eos().template read<bool>(); }
  
  private:
    struct
    {
      binary_p eoc;
    } params;
    array<stage> stages;
    size_t stageIdx;
    stage final;
  };

  template<typename T>
  class ad : public envelope<T>
  {
    typename envelope<T>::stage attackStage;

  public:
    ad(analog_t attackDuration, analog_t decayDuration, analog_t sampleRate) : envelope<T>(&attackStage, 1, sampleRate), attackStage(sampleRate)
    { attack().target() = T(1); attack().duration() = attackDuration; decay().duration() = decayDuration; }
    
    typename envelope<T>::stage& attack() { return attackStage; }
    typename envelope<T>::stage& decay() { return envelope<T>::finalStage(); }
    using envelope<T>::eoc;

    using envelope<T>::trigger;
    using envelope<T>::generate;
    using envelope<T>::elementAt;
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
    void trigger() override { attack().target() = T(1); ad<T>::trigger(); gateOn = false; }
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
      attackStage.duration() = attackDuration;
      attackStage.target() = T(1);
      decayStage.duration() = decayDuration;
      decayStage.target() = sustainLevel;
      envelope<T>::finalStage().duration() = releaseDuration;
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
  class slew : public unitProcessor<T>, protected plist<5>
  {
  public:
    // note: choice of epsilon will depend on the amount of noise in the signal to be slewed.
    // the default value was chosen based on testing with an OWL module's audio input.
    slew(analog_t sampleRate, analog_t riseRate, analog_t fallRate, T initialValue = T(0), T epsilon = math::epsilon<T>()*1000)
    : unitProcessor<T>(), epsilon(epsilon), dt(1.0f/sampleRate)
    {
      params.rise.value = riseRate; params.fall.value = fallRate; params.output.value = initialValue;
    }
    
    void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate; }
    const parameters& getParameters() const override { return *this; }

    parameter rise() const { return params.rise({ "rise", 'a', analog_p::type }); }
    parameter fall() const { return params.fall({ "fall", 'd', analog_p::type }); }
    parameter rising() const { return params.rising({ "rising", 'r', binary_p::type }); }
    parameter falling() const { return params.falling({ "falling", 'f', binary_p::type }); }
    parameter value() const { return params.output({ "value", 'v', param<T>::type }); }
    
    T process(const T& v) override;
    using processor<T>::process;
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter plist[plsz] = { rise(), fall(), rising(), falling(), value() };
      return plist[index];
    }
    
  private:
    struct
    {
      analog_p rise;
      analog_p fall;
      binary_p rising;
      binary_p falling;
      param<T> output;
    } params;
    T epsilon;
    analog_t dt;
  };

  // note: W must implement waveform<T>
  template<class W>
  class oscil final : public unitGenerator<typename W::SampleType>, protected plist<4>
  {
  public:
    using T = typename W::SampleType;
    
    oscil() : unitGenerator<T>(), phase(0), dt(0) { params.fHz.value = 440.0; }

    template<typename... Ts>
    explicit oscil(analog_t sampleRate, analog_t freqInHz, Ts... wargs)
    : unitGenerator<T>(), waveform(wargs...), phase(0), dt(1.0f/sampleRate)
    { params.fHz.value = freqInHz; }
    
    void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate;}

    W waveform;
    
    const parameters& getParameters() const override { return *this; }
    
    // frequency in Hz without FM applied
    parameter fHz() const { return params.fHz({ "frequency", 'f', analog_p::type }); }
    // linear frequency modulation
    parameter fmLin() const { return params.fmLin({ "fm (lin)", 'l', analog_p::type }); }
    // v/oct (exponential) frequency modulation
    parameter fmExp() const { return params.fmExp({ "fm (v/oct)", 'v', analog_p::type }); }
    // phase modulation
    parameter pm() const { return params.pm({ "phase mod", 'p', analog_p::type }); }

    T generate() override;

    void reset() { phase = 0; }
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] { fHz(), fmLin(), fmExp(), pm() };
      return p[index];
    }
    
  private:
    struct
    {
      analog_p fHz;
      analog_p fmLin;
      analog_p fmExp;
      analog_p pm;
    } params;
    analog_t phase;
    analog_t dt;
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
  class delay : public unitProcessor<T>, protected plist<2>
  {
  public:
    delay(array<T> buffer, analog_t sampleRate, analog_t delayInSeconds = 0, analog_t feedbackAmount = 0)
    : unitProcessor<T>(), buffer(buffer.getData(), buffer.getSize()), dt(1.0f/sampleRate)
    {
      params.time.value = duration::fromSeconds(delayInSeconds, sampleRate);
      mDelayInSamples = params.time.value.samples;
      params.feedback.value = feedbackAmount;
    }
    
    void setSampleRate(float sampleRate) override { dt = 1.0f / sampleRate; }
    const parameters& getParameters() const override { return *this; }

    delayline<T>& getBuffer() { return buffer; }
    const delayline<T>& getBuffer() const { return buffer; }

    /// delay time expressed as vessl::duration (i.e. samples), can be set using an analog_t
    parameter time() const { return params.time({ "time", 't', param<duration>::type }); }
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    parameter feedback() const { return params.feedback({ "feedback", 'f', analog_p::type }); }
    
    T process(const T& in) override;

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<duration::mode TimeMode = duration::mode::slew>
    void process(array<T> input, array<T> output);
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] = { time(), feedback() }; 
      return p[index];
    }

  private:
    struct
    {
      param<duration> time;
      analog_p feedback;
    } params;
    
  protected:
    delayline<T> buffer;
    analog_t mDelayInSamples;
    analog_t dt;
  };

  template<typename T>
  class follow : public unitProcessor<T>, protected plist<1>
  {
  public:
    follow(array<T> window, analog_t sampleRate, analog_t responseTimeInSeconds)
    : unitProcessor<T>(), writer(window), window(window)
    , delta(math::exp(-1.0 / (sampleRate*responseTimeInSeconds)))
    , previous(0), current(0)
    {
      params.response.value = responseTimeInSeconds;
    }
    
    void setSampleRate(float sampleRate) override { delta = math::exp(-1.0 / (sampleRate*params.response.value)); }
    const parameters& getParameters() const override { return *this; }
    
    parameter response() const { return params.response({ "response time", 'r', analog_p::type }); }

    T process(const T& in) override;

    using processor<T>::process;
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] = { response() }; return p[index];
    }
    
  private:
    struct
    {
      // @todo actually use this parameter
      analog_p response;
    } params;
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
  class freeze : public unit, public processor<T>, public generator<T>, protected plist<4>
  {
  public:
    explicit freeze(array<T> buffer, float sampleRate) : unit(), buffer(buffer.getData(), buffer.getSize())
    , mPhase(0), crossfade(0.75f), freezeDelay(0), freezeSize(0), readRate(1), dt(1.0/sampleRate)
    {
      params.size.value.samples = buffer.getSize()-1;
      freezeSize = params.size.value.samples;
      params.rate.value = 1.0;
    }
    
    void setSampleRate(analog_t sr) override { dt = 1.0f/sr; }
    const parameters& getParameters() const override { return *this; }
    
    delayline<T>& getBuffer() { return buffer; }
    const delayline<T>& getBuffer() const { return buffer; }
    
    parameter enabled() const { return params.enabled({ "enabled", 'e', binary_p::type }); }
    // end of the freeze loop in samples relative to the most recently recorded sample
    parameter position() const { return params.position({ "position", 'p', analog_p::type }); }
    // size of the freeze loop as a duration (samples).
    // the beginning of the freeze loop, when played forward, will be position + size.
    parameter size() const { return params.size({ "size", 's', analog_p::type }); }
    // rate of playback when enabled, can be negative to play in reverse
    parameter rate() const { return params.rate({ "rate", 'r', analog_p::type }); }
    // should this be a parameter? there's not much gained by it.
    analog_t phase() const { return mPhase; }
    // reset the phase to zero (argument for making it a parameter?)
    void reset() { mPhase = 0; }

    T generate() override;

    T process(const T& in) override;

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<duration::mode TimeMode = duration::mode::slew>
    void generate(array<T> output) { procgen<TimeMode, false>(output, output); }

    template<duration::mode TimeMode = duration::mode::slew>
    void process(array<T> input, array<T> output) { procgen<TimeMode, true>(input, output); }
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] = { enabled(), position(), size(), rate() };
      return p[index];
    }

  private:
    // shared routine used by templated processing and generation methods
    template<duration::mode TimeMode, bool UseInput>
    void procgen(array<T> input, array<T> output);
    
    struct
    {
      binary_p enabled;
      analog_p position;
      param<duration> size;
      analog_p rate;
    } params;
    delayline<T> buffer;
    analog_t mPhase;
    smoother<analog_t> crossfade; // used to crossfade between the incoming signal and the freeze signal when enabled changes.
    analog_t freezeDelay, freezeSize, readRate, dt;
  };

  // unit for use with filter types defined in the filtering namespace.
  // For example, for a 2 stage biquad low pass filter use:
  // filter<float, filtering::biquad<2>::lowPass>
  template<typename T, template<typename> typename H>
  class filter : public unitProcessor<T>, protected plist<3>
  {
  public:
    typedef H<T> function;
    
    explicit filter(analog_t sampleRate) : unitProcessor<T>(), sampleRate(sampleRate)
    {
      params.fHz.value = sampleRate;
      params.q.value = 1;
      params.emphasis.value = gain::fromDecibels(0);
    }
    
    filter(analog_t sampleRate, analog_t freqInHz, analog_t kyu = filtering::q::butterworth<analog_t>(), gain emphasis = gain::fromDecibels(0) ) 
    : unitProcessor<T>(), sampleRate(sampleRate)
    {
      params.fHz.value = freqInHz;
      params.q.value = kyu;
      params.emphasis.value = emphasis;
    }
    
    void setSampleRate(analog_t sampleRate) override { this->sampleRate = sampleRate; }
    const parameters& getParameters() const override { return *this; }

    parameter fHz() const { return params.fHz({ "fHz", 'f', analog_p::type }); }
    parameter q() const { return params.q({ "q", 'q', analog_p::type }); }
    // unused by some filter types (see filtering section)
    parameter emphasis() const { return params.emphasis({ "emphasis", 'e', analog_p::type }); }
    
    T process(const T& in) override
    {
      T out;
      filtering::args fargs = {
        sampleRate,
        params.fHz.value,
        math::max(params.q.value, 0.01),
        params.emphasis.value
      };
      func.process(&in, &out, 1, fargs);
      return out;
    }
    
    void process(array<T> in, array<T> out) override
    {
      filtering::args fargs = {
        sampleRate,
        params.fHz.value,
        math::max(params.q.value, 0.01),
        params.emphasis.value
      };
      func.process(in.getData(), out.getData(), in.getSize(), fargs);
    }
    
  protected:
    parameter elementAt(size_t index) const override { parameter p[plsz] = {fHz(), q(), emphasis()}; return p[index]; }
    
  private:
    struct
    {
      analog_p fHz;
      analog_p q;
      param<gain> emphasis;
    } params;
    function func;
    analog_t sampleRate;
  };

  // designed to work with floating point types.
  template<typename T, uint32_t MaxBits>
  class bitcrush : public unitProcessor<T>, protected plist<3>
  {
  public:
    bitcrush(analog_t sampleRate, analog_t bitRate, analog_t bitDepth = MaxBits)
      : unitProcessor<T>(), prevInput(0), currSample(0), rateAlpha(0), dt(1.0f/sampleRate)
    {
      params.bitRate.value = bitRate;
      params.bitDepth.value = bitDepth;
    }
    
    void setSampleRate(analog_t sampleRate) override { dt = 1.0f / sampleRate;}
    const parameters& getParameters() const override { return *this; }

    parameter rate() const { return params.bitRate({ "bit rate", 'r', analog_p::type }); }
    parameter depth() const { return params.bitDepth({ "bit depth", 'd', analog_p::type }); }
    parameter mangle() const { return params.mangle({ "mangle", 'm', binary_p::type }); }

    T process(const T& in) override;

    using unitProcessor<T>::process;
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] = { rate(), depth(), mangle() };
      return p[index];
    }
    
  private:
    struct
    {
      analog_p bitRate;
      analog_p bitDepth;
      binary_p mangle;
    } params;
    T prevInput;
    T currSample;
    analog_t rateAlpha;
    analog_t dt;
  };

  // A simple peak limiter adapted from pinchenettes/stmlib via DaisySP
  template<typename T>
  class limiter : public unitProcessor<T>, protected plist<2>
  {
  public:
    explicit limiter(gain preGain = gain::fromDecibels(0)) : unitProcessor<T>()
    {
      params.preGain.value = preGain;
      params.peak.value = 0.5;
    }
    
    const parameters& getParameters() const override { return *this; }

    parameter preGain() const { return params.preGain({ "pre-gain", 'g', param<gain>::type }); }
    parameter peak() const { return params.peak({ "peak", 'k', param<T>::type }); }

    T process(const T& in) override;
    using unitProcessor<T>::process;
    
  protected:
    parameter elementAt(size_t index) const override
    {
      parameter p[plsz] = { preGain(), peak() };
      return p[index];
    }
    
  private:
    struct
    {
      param<gain> preGain;
      param<T> peak;
    } params;
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

  // @todo ARM specialization
  template<typename T>
  void array<T>::writer::write(const reader& r)
  {
    size_t rsz = r.available();
    size_t wsz = available();
    VASSERT(wsz >= rsz, "Not enough space in writer for the contents of reader");
    const T* rh = *r;
    memcpy(static_cast<void*>(head), static_cast<const void*>(rh), rsz * sizeof(T));
    head += rsz;
  }
  
  template<typename T>
  void array<T>::copyTo(array dest)
  {
    writer w(dest);
    reader r(*this);
    w.write(r);
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
  array<T> array<T>::offset(T value, array dest) const
  {
    VASSERT(size <= dest.size, "arrays are have different lengths or destination is too small");
    reader a(data, size);
    writer b(dest.data, dest.size);
    while (a)
    {
      b << a.read() + value;
    }
    return dest;
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
  
  template<size_t STAGES>
  template<typename T, class CoGen>
  void filtering::biquad<STAGES>::df2T<T, CoGen>::process(const T* source, T* dest, size_t blockSize, const args& args)
  {
    // update our coefficents
    cg(coeff.getData(), args.omega(), args.q, args.g);
    
    // run the filter
    const T* input = source;
    array<analog_t>::reader cor = coeff.getReader();
    for (size_t s = 0; s < STAGES; s++)
    {
      analog_t b0 = cor.read(); analog_t b1 = cor.read(); analog_t b2 = cor.read();
      analog_t a1 = cor.read();  analog_t a2 = cor.read();
      T* st = state.getData() + 2*s;
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

  template<size_t STAGES>
  template<typename T>
  void filtering::biquad<STAGES>::copy(T* coeff) 
  {
    if (STAGES > 1)
    {
      array<T> src(coeff, COEFF_NUM);
      for (size_t i = 1; i < STAGES; i++)
      {
        array<T> dst(coeff + COEFF_NUM*i, COEFF_NUM);
        src.copyTo(dst);
      }
    }
  }

  template<size_t STAGES>
  void filtering::biquad<STAGES>::lpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::hpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::bpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::ntcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain _) const
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::pkcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.toDecibels())/20);
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::lscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.toDecibels())/20);
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

  template<size_t STAGES>
  void filtering::biquad<STAGES>::hscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.toDecibels())/20);
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

  template<typename T, typename N>
  T noiseGenerator<T, N>::generate()
  {
    step += dt * rate();
    if (step >= 1)
    {
      value = next;
      next = noise();
      step = math::wrap01(step);
    }
    return value;
  }

  template<typename T, typename N>
  template<typename E>
  T noiseGenerator<T, N>::generate()
  {
    T s = generate();
    return easing::interp<E>(s, next, step);
  }

  template<typename T>
  void ramp<T>::trigger()
  {
    if (duration() > math::epsilon<analog_t>())
    {
      t = 0;
      params.eor.value = false;
    }
    else
    {
      t = 1;
      params.eor.value = true;
    }
  }

  template<typename T>
  T ramp<T>::generate()
  {
    analog_t lt = t;
    if (isActive())
    {
      analog_t dinv = 1.f / duration();
      t += dt * dinv;
      if (t >= 1)
      {
        params.eor.value = true;
      }
    }
    return easing::lerp(params.from.value, params.to.value, lt);
  }

  template<typename T>
  template<typename E>
  T envelope<T>::stage::step()
  {
    analog_t s = dt / math::max<analog_t>(static_cast<analog_t>(duration()), dt);
    time += s;
    analog_t t = math::constrain<analog_t>(time, 0.0, 1.0);
    params.output.value = easing::interp<E, T>(begin, params.target.value, t);
    if (time >= 1)
    {
      params.active.value = false;
      params.eos.value = true;
    }
    return params.output.value;
  }

  template<typename T>
  void envelope<T>::trigger()
  {
    for (stage& stage : stages)
    {
      stage.reset();
    }
    final.reset();
    params.eoc.value = false;
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
      params.eoc.value = true;
    }
    else if (shouldAdvance(stageIdx))
    {
      getStage(++stageIdx).start(value);
    }
    return value;
  }

  template<typename T>
  void envelope<T>::setSampleRate(analog_t sampleRate)
  {
    for (stage& stage : stages)
    {
      stage.setSampleRate(sampleRate);
    }
    final.setSampleRate(sampleRate);
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
    attack().target() = math::max(value, attackTarget);
  }

  template<typename T>
  T slew<T>::process(const T& v)
  {
    binary_t isRise = v > params.output.value+epsilon;
    binary_t isFall = v < params.output.value-epsilon;
    if (isRise)
    {
      params.output.value = math::min(v, params.output.value + params.rise.value*dt);
    }
    else if (isFall)
    {
      params.output.value = math::max(v, params.output.value - params.fall.value*dt);
    }
    else
    {
      params.output.value = v;
    }
    params.rising.value = isRise;
    params.falling.value = isFall;
    return params.output.value;
  }

  template<class W>
  typename W::SampleType oscil<W>::generate()
  {
    typename W::SampleType val = waveform.evaluate(phase + pm());
    phase += (fHz() * math::exp2(static_cast<analog_t>(fmExp())) + fmLin()) * dt;
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
    mDelayInSamples = params.time.value.samples;
    // delay time in samples
    analog_t dts = math::constrain<analog_t>(mDelayInSamples, 0.f, static_cast<analog_t>(buffer.getSize()-1));
    analog_t s = buffer.template read<I>(dts);
    analog_t fbk = math::constrain<analog_t>(static_cast<analog_t>(feedback()), -1.0, 1.0);
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
      analog_t dst = dt * 10.0f;
      while (r && w)
      {
        T in = r.read();
        mDelayInSamples = easing::lerp(mDelayInSamples, params.time.value.samples, dst);
        // delay time in samples
        analog_t dts = math::constrain<analog_t>(mDelayInSamples, 0.f, static_cast<analog_t>(buffer.getSize()-1));
        analog_t wet = buffer.template read<I>(dts);
        analog_t fbk = math::constrain<analog_t>(static_cast<analog_t>(feedback()), -1.0, 1.0);
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
      analog_t targetSampleDelay = params.time.value.samples;
      // delay time in samples
      analog_t fts = math::constrain<analog_t>(mDelayInSamples, 0.0, static_cast<analog_t>(buffer.getSize()-1));
      analog_t tts = math::constrain<analog_t>(targetSampleDelay, 0.0, static_cast<analog_t>(buffer.getSize()-1));
      analog_t fbk = math::constrain<analog_t>(static_cast<analog_t>(feedback()), -1.0, 1.0);
      
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

  template<typename T>
  T follow<T>::process(const T& in) 
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

  template<typename T, typename I>
  T freeze<T, I>::generate() 
  {
    freezeDelay = static_cast<analog_t>(position());
    freezeSize  = params.size.value.samples;
    analog_t sampleDelay = freezeDelay + (1.0-mPhase)*freezeSize;
    mPhase = math::wrap01(mPhase + rate() / freezeSize);
    return buffer.template read<I>(sampleDelay);
  }

  template<typename T, typename I>
  T freeze<T, I>::process(const T& in) 
  {
    binary_t isEnabled = static_cast<binary_t>(enabled());
    analog_t wetLevel = crossfade = (isEnabled ? 1.0 : 0.0);
    T wet = generate();
    if (!isEnabled)
    {
      buffer.write(in);
    }
    return wetLevel*wet + (1.0 - wetLevel)*in;
  }

  template<typename T, typename I>
  template<duration::mode TimeMode, bool UseInput>
  void freeze<T, I>::procgen(array<T> input, array<T> output) 
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
      analog_t fd = params.position.value;
      analog_t fs = params.size.value.samples;
      analog_t rt = params.rate.value;
      analog_t st = dt;
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
          binary_t isEnabled = params.enabled.value;
          analog_t wetLevel = crossfade = (isEnabled ? 1.0 : 0.0);
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
      analog_t fd0 = freezeDelay,   fd1 = params.position.value;
      analog_t fs0 = freezeSize,    fs1 = params.size.value.samples;
      analog_t fade = 0, fadeInc = 1.0 / output.getSize();
      analog_t r0 = readRate, r1 = params.rate.value;
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
          binary_t isEnabled = params.enabled.value;
          analog_t wetLevel = crossfade = (isEnabled ? 1.0 : 0.0);
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

  template<typename T, uint32_t MaxBits>
  T bitcrush<T, MaxBits>::process(const T& in)
  {
    rateAlpha += math::max<analog_t>(1.0, static_cast<analog_t>(rate()))*dt;
    if (rateAlpha >= 1)
    {
      rateAlpha -= 1;
      currSample = easing::lerp(prevInput, in, rateAlpha);
    }

    analog_t bd = math::constrain<analog_t>(static_cast<analog_t>(depth()), 2.0, MaxBits);
    analog_t scalar = math::pow<analog_t>(2.0, bd) - 1;
    T val = math::round(currSample*scalar);
    if (mangle())
    {
      val = math::xore(val, math::round(prevInput*scalar));
    }
    prevInput = in;
    return val * (1.0 / scalar);
  }

  template<typename T>
  T limiter<T>::process(const T& in)
  {
    T pre = in*params.preGain.value.toScale();
    T peak = math::abs(pre);
    T error = peak - params.peak.value;
    params.peak.value += (error > T(0)) ? T(0.05) : T(0.00002) * error;
    analog_t gain = (params.peak.value <= 1.0 ? 1.0 : 1.0 / params.peak.value);
    // DaisySP returns this, which sounds better for how I typically use this.
    return saturation::softlimit(pre*gain*0.7);
    // stmlib returns this, which clips more easily, but is faster.
    //return pre*gain*0.8;
  }

#ifdef ARM_CORTEX
  template<>
  void array<float>::fill(float value)
  {
    arm_fill_f32(value, data, size);  
  }

  template<>
  void array<float>::copyTo(array dest)
  {
    VASSERT(this->getSize() <= dest.getSize(), "Not enough room in destination for this array");
    arm_copy_f32(data, dest.getData(), this->getSize());
  }
  
  template<>
  array<float> array<float>::offset(float value, array dest) const
  {
    VASSERT(this->getSize() <= dest.getSize(), "Not enough room in destination for this array");
    arm_offset_f32(data, value, dest.getData(), this->getSize());
    return dest;
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
  
  namespace math
  {
    template<>
    inline float sin(float x) { return arm_sin_f32(x); }

    template<>
    inline float sqrt(float x)
    {
      float out;
      if (ARM_MATH_SUCCESS == arm_sqrt_f32(x, &out))
      {
        return out;
      }
      return 0;
    }
  }
  
  namespace filtering
  {
    template<size_t STAGES>
    template<class CoGen>
    struct biquad<STAGES>::df2T<float, CoGen> final : cascade<float, 2>
    {
      arm_biquad_cascade_df2T_instance_f32 inst;
          
      using biquad<STAGES>::cascade<float, 2>::coeff;
      using biquad<STAGES>::cascade<float, 2>::state;
      using biquad<STAGES>::cascade<float, 2>::getCoeffSize;
      using biquad<STAGES>::cascade<float, 2>::getStateSize;
    
      static CoGen cg;
  
      df2T() { arm_biquad_cascade_df2T_init_f32(&inst, STAGES, coeff.getData(), state.getData()); }
    
      size_t getStageCount() const { return STAGES; }
  
      void process(const float* source, float* dest, size_t blockSize, const args& args)
      {
        cg(coeff.getData(), args.omega(), args.q, args.g);
        arm_biquad_cascade_df2T_f32(&inst, source, dest, blockSize);
      }
    };
  }
#endif
}
