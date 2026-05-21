////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2025-2026 Damien Quartz
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
#include <cstring>
#include <limits>
#include <utility>

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

#ifdef cos
#undef cos
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

// MSVC check: Source - https://stackoverflow.com/a/77012222
// Posted by Peter Cordes
// Retrieved 2026-05-17, License - CC BY-SA 4.0
#if defined(_MSC_VER) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#define VESSL_INLINE __forceinline
#else
#define VESSL_INLINE __attribute__((always_inline)) inline
#endif


// Note: In all classes using typename T, it is assumed to be POD and to have support for all arithmetic operators
namespace vessl
{
  using char_t    = char;
  using size_t    = uint64_t;
  using binary_t  = bool;
  using digital_t = int64_t;
  using analog_t  = float;
  using phase_t   = uint32_t;

  static constexpr phase_t phase_360  = UINT32_MAX;
  static constexpr phase_t phase_180  = UINT32_MAX >> 1;
  static constexpr phase_t phase_90   = phase_180 >> 1;
  static constexpr phase_t phase_270  = phase_180 + phase_90;
  static constexpr phase_t phase_zero = 0UL;

  // we use this in place of static_cast throughout the library for non-pointer types
  // so that we can specialize conversions between some of our value types (e.g. phase_t <--> analog_t)
  template<typename T, typename F>
  VESSL_INLINE constexpr T cast(const F& from) { return static_cast<T>(from); }

  // phase t <--> digital_t (@todo convert to degrees?)
  template<>
  VESSL_INLINE constexpr digital_t cast<digital_t, phase_t>(const phase_t& from)
  {
    return from > phase_180 ? 1 : 0;
  }

  template<>
  VESSL_INLINE constexpr phase_t cast<phase_t, digital_t>(const digital_t& from)
  {
    return from <= 0 ? phase_zero : phase_360;
  }

  template<>
  VESSL_INLINE constexpr analog_t cast<analog_t, phase_t>(const phase_t& from)
  {
    return static_cast<analog_t>(from) / 4294967295.0f;
  }

  template<>
  VESSL_INLINE constexpr phase_t cast<phase_t, analog_t>(const analog_t& from)
  {
    return static_cast<phase_t>(from * 4294967295.0);
  }

  template<>
  VESSL_INLINE constexpr phase_t cast<phase_t, size_t>(const size_t& from)
  {
    return from % phase_360;
  }
  
  // phase_t <--> binary_t
  template<>
  VESSL_INLINE constexpr binary_t cast<binary_t, phase_t>(const phase_t& from) { return from > phase_180; }
  
  template<>
  VESSL_INLINE constexpr phase_t cast<phase_t, binary_t>(const binary_t& from) { return from ? phase_360 : phase_zero; }
  
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

    [[nodiscard]] virtual binary_t is_empty() const = 0;
    virtual T read() = 0;

    VESSL_INLINE explicit operator bool() const { return !is_empty(); }
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

    virtual binary_t is_full() const = 0;
    virtual void write(const T& value) = 0;

    VESSL_INLINE sink& operator<<(const T& value) { write(value); return *this; }
    VESSL_INLINE explicit operator binary_t() const { return !is_full(); }
  };

  template<typename T>
  class generator : public source<T>
  {
  public:
    generator() = default;
    virtual T generate() = 0;

    // by default, we assume an endless source of data
    VESSL_INLINE binary_t is_empty() const override { return false; }
    VESSL_INLINE T read() override { return generate(); }
  };
  
  template<typename T>
  class list
  {
  public:
    virtual ~list() = default;
    virtual size_t size() const = 0;
    VESSL_INLINE binary_t is_empty() const { return size() == 0; }
    VESSL_INLINE T operator[](size_t index) const
    {
      VASSERT(index < size(), "Attempt to access a list element with out-of-bounds index."); 
      return element_at(index);
    }
    
    class iterator
    {
    public:
      static iterator begin(const list& src) { return iterator(src, 0); }
      static iterator end(const list& src) { return iterator(src, src.size()); }
      T operator*() const { return (*src_)[index_]; }
      iterator& operator++() { index_++; return *this; } 
      bool operator==(const iterator& it) const { return src_ == it.src_ && index_ == it.index_; }
      bool operator!=(const iterator& it) const { return !(*this == it); }
      
    private:
      const list* src_ = nullptr;
      size_t index_ = 0;
      iterator(const list& source, size_t index): src_(&source), index_(index) {}
    };
    
  protected:
    virtual T element_at(size_t index) const = 0;
  };
  
  template<typename T>
  typename list<T>::iterator begin(const list<T>& lst) { return list<T>::iterator::begin(lst); }
  
  template<typename T>
  typename list<T>::iterator end(const list<T>& lst) { return list<T>::iterator::end(lst); }

  template<typename T>
  class array // not a list<T> to keep size to 16 bytes
  {
  public:
    VESSL_INLINE array() : data_(nullptr), size_(0) {}
    VESSL_INLINE array(T* src_data, size_t src_size) : data_(src_data), size_(src_size) {}

    VESSL_INLINE T* data() { return data_; }
    VESSL_INLINE const T* data() const { return data_; }
    VESSL_INLINE size_t size() const { return size_; }

    // ranged-based for support
    VESSL_INLINE T* begin() { return data_; }
    VESSL_INLINE T* end() { return data_ + size_; }
    
    VESSL_INLINE T& operator[](size_t index) { return data_[index]; }
    VESSL_INLINE const T& operator[](size_t index) const { return data_[index]; }

    class reader final : public source<T>
    {
    public:
      VESSL_INLINE explicit reader() : source<T>(), begin_(nullptr), head_(nullptr), end_(nullptr) {}
      VESSL_INLINE explicit reader(array source) : source<T>(), begin_(source.data_), head_(source.data_), end_(source.data_ + source.size_) {}
      VESSL_INLINE reader(const T* data, size_t size) : source<T>(), begin_(data), head_(data), end_(data + size) {}

      // source methods
      VESSL_INLINE binary_t is_empty() const override { return head_ == end_; }
      VESSL_INLINE T read() override { return *head_++; }

      VESSL_INLINE size_t available() const { return end_ - head_; }

      VESSL_INLINE T peek() const { return *head_; }
      VESSL_INLINE const T* operator*() const { return head_; }

      VESSL_INLINE reader reset() { head_ = begin_; return *this; }
      
    protected:
      const T* begin_;
      const T* head_;
      const T* end_;
    };

    class writer final : public sink<T>
    {
    public:
      VESSL_INLINE explicit writer(array source) : sink<T>(), head_(source.data_), end_(source.data_ + source.size_) {}
      VESSL_INLINE writer(T* data, size_t size) : sink<T>(), head_(data), end_(data + size) {}

      VESSL_INLINE binary_t is_full() const override { return head_ == end_; }
      VESSL_INLINE void write(const T& v) override { *head_++ = v; }
      // block copy the entire contents of reader into this writer.
      // writer must have enough space for the contents of reader.
      // ReSharper disable once CppEnforceOverridingFunctionStyle
      void write(const reader& r);
      VESSL_INLINE size_t available() const { return end_ - head_; }
      
    protected:
      T* head_;
      const T* end_;
    };

    [[nodiscard]] VESSL_INLINE reader make_reader() const { return reader(*this); }
    [[nodiscard]] VESSL_INLINE writer make_writer() { return writer(*this); }

    // block copy this array to dest, which must be large enough to hold this array.
    void copy_to(array dest);

    void fill(T value);
    
    // adds value to every element in the array, returns dest
    array offset(T value, array dest) const;
    VESSL_INLINE array offset(T value) { return offset(value, *this); }
    
    // element-wise addition of this and other, returns dest
    array add(array other, array dest) const;
    VESSL_INLINE array add(array other) { return add(other, *this); }
    
    // element-wise subtraction of this and other, returns dest
    array subtract(array other, array dest) const;
    VESSL_INLINE array subtract(array other) { return subtract(other, *this); }
    
    // scales every element in this array by value, returns dest
    array scale(T value, array dest) const;
    VESSL_INLINE array scale(T value) { return scale(value, *this); }
    
    // element-wise multiplication of this and other, returns dest
    array multiply(array other, array dest) const;
    VESSL_INLINE array multiply(array other) { return multiply(other, *this); }
    
  protected:
    T* data_;
    size_t size_;
  };

  template<typename T>
  VESSL_INLINE T* begin(array<T>& arr) { return arr.begin(); }

  template<typename T>
  VESSL_INLINE T* end(array<T>& arr) { return arr.end(); }
  
  // for the most part these resolve to standard math calls,
  // template specializations are provided for ARM_CORTEX where applicable.
  namespace math
  {
    template<typename T>
    VESSL_INLINE constexpr T e() { return cast<T>(2.71828182845904523536); }
    
    template<typename T>
    VESSL_INLINE constexpr T pi() { return cast<T>(3.1415926535897932385); }

    template<>
    VESSL_INLINE constexpr phase_t pi<phase_t>() { return phase_180; }

    template<typename T>
    VESSL_INLINE constexpr T two_pi() { return pi<T>() * 2; }

    template<>
    VESSL_INLINE constexpr phase_t two_pi() { return phase_360; }

    template<typename T>
    VESSL_INLINE T abs(const T& val) { return ::abs(val); }
    
    template<typename T>
    VESSL_INLINE T constrain(T val, T low, T high) { return val < low ? low : val > high ? high : val; }
    
    template<typename T>
    VESSL_INLINE T epsilon() { return std::numeric_limits<T>::epsilon(); }

    template<typename T>
    VESSL_INLINE T exp(T v) { return ::exp(v); }

    template<typename T>
    VESSL_INLINE T exp2(T v) { return ::exp2(v); }

    template<typename T>
    VESSL_INLINE T exp10(T v) { return ::pow(T(10), v); }
    
    template<typename T>
    VESSL_INLINE T log(T v) { return ::log(v); }

    template<typename T>
    VESSL_INLINE T log10(T v) { return ::log10(v); }
    
    template<typename T>
    VESSL_INLINE T max(const T& a, const T& b) { return a > b ? a : b; }

    template<typename T>
    VESSL_INLINE T min(const T& a, const T& b) { return a < b ? a : b; }

    template<typename T>
    VESSL_INLINE T mod(T v, T* i) { return std::modf(v, i); }

    template<typename T>
    VESSL_INLINE T pow(T x, T y) { return ::pow(x, y); }

    template<typename T>
    VESSL_INLINE T round(T x) { return ::round(x); }
    
    template<typename T>
    VESSL_INLINE T floor(T x) { return ::floor(x); }

    template<typename T, typename R = T>
    VESSL_INLINE T sin(R r) { return std::sin(r); }

    template<typename T, typename R = T>
    VESSL_INLINE T cos(R r) { return std::cos(r); }

    template<typename T>
    VESSL_INLINE T sqrt(T x) { return std::sqrt(x); }
    
    template<typename T>
    VESSL_INLINE T sqrt2() { static T v = sqrt(2); return v; }

    template<typename T>
    VESSL_INLINE T tan(T x) { return std::tan(x); }
    
    template<typename T>
    VESSL_INLINE T wrap(T val, T low, T high)
    {
      // @todo probably a way to do this without while loops.
      T diff = high - low;
      while (val < low) { val += diff; }
      while (val > high) { val -= diff; }
      return val;
    }

    // @todo use modf here
    template<typename T>
    VESSL_INLINE T wrap01(T val) { return wrap(val, T(0), T(1)); }

    template<>
    VESSL_INLINE phase_t wrap01<phase_t>(phase_t val) { return val; }
    
    template<>
    VESSL_INLINE analog_t wrap01<analog_t>(analog_t v) { analog_t i; analog_t f = mod(v, &i); return f < 0 ? f + 1.0f : f; }

    template<typename T>
    VESSL_INLINE binary_t is_nan(T n) { return isnan(n); }
    
    // xore because xor is a keyword
    template<typename T>
    VESSL_INLINE T xore(const T& a, const T& b)
    {
      return a ^ b;
    }
    
    template<>
    VESSL_INLINE analog_t xore(const analog_t& a, const analog_t& b)
    {
      return xore(cast<digital_t>(a), cast<digital_t>(b));
    }
  
    template<typename T>
    struct matrix_data
    {
      T*       data;
      uint32_t num_rows;
      uint32_t num_cols;
      
      VESSL_INLINE matrix_data() : data(nullptr), num_rows(0), num_cols(0) {}
      VESSL_INLINE matrix_data(T* d, uint32_t r, uint32_t c) : data(d), num_rows(r), num_cols(c) {}
      
      [[nodiscard]] VESSL_INLINE T* operator*() { return data; }
      [[nodiscard]] VESSL_INLINE T* operator*() const { return data; }
      [[nodiscard]] VESSL_INLINE uint32_t rows() const { return num_rows; }
      [[nodiscard]] VESSL_INLINE uint32_t cols() const { return num_cols; }
    };
    
    template<typename T>
    class matrix
    {
    public:
      VESSL_INLINE matrix() = default;
      VESSL_INLINE matrix(T* src_data, size_t rows, size_t cols) : data_(src_data, rows, cols) {}

      VESSL_INLINE T* data() { return *data_; }
      VESSL_INLINE const T* data() const { return *data_; }
      [[nodiscard]] VESSL_INLINE size_t rows() const { return data_.rows(); }
      [[nodiscard]] VESSL_INLINE size_t columns() const { return data_.cols(); }
      [[nodiscard]] VESSL_INLINE size_t size() const { return rows()*columns(); }
      
      VESSL_INLINE T* operator[](uint32_t row) { return &data()[row*columns()]; }
      VESSL_INLINE const T* operator[](uint32_t row) const { return &data()[row*columns()]; }
      
      VESSL_INLINE void clear() { array<T> arr(data(), size()); arr.fill(T(0LL)); }
      VESSL_INLINE T get(size_t row, size_t col) const { return data()[row*columns() + col]; }
      VESSL_INLINE void set(size_t row, size_t col, T value) { data()[row*columns() + col] = value; }
      
      // element-wise addition of this and other, returns dest
      matrix add(matrix other, matrix dest) const;
      VESSL_INLINE matrix add(matrix other) { return add(other, *this); }
      
      // element-wise subtraction of this and other, returns dest
      matrix subtract(matrix other, matrix dest) const;
      VESSL_INLINE matrix subtract(matrix other) { return subtract(other, *this); }
      
      // scales every element in this matrix by value, returns dest
      matrix scale(T value, matrix dest) const;
      VESSL_INLINE matrix scale(T value) { return scale(value, *this); }
      
      // matrix multiplication of this and other, returns dest
      matrix multiply(matrix other, matrix dest) const;
      VESSL_INLINE matrix multiply(matrix other) { return multiply(other, *this); }
      
      // matrix vector multiplication of this and vector, returns dest
      array<T> multiply(const array<T>& vector, array<T> dest) const;
      
    private:
      matrix_data<T> data_;
    };
  }

  namespace time
  {
    // can be used by units as an indication for how to treat changes in time values (see: delay & freeze)
    enum class mode : uint8_t
    {
      snap, // use duration value directly
      slew, // smooth duration when it changes
      fade, // crossfade to new duration across a block of samples
    };
  }

  // @todo move built-in types like analog_t into this namespace?
  namespace types
  {
    // stored as decibels (0 = unity gain).
    struct gain
    {
      analog_t db;
      
      gain() : db(0) {}
      explicit gain(analog_t dbv) : db(dbv) {}
      static analog_t decibels_to_scale(analog_t db) { return math::exp10(db*cast<analog_t>(0.05));}
      static analog_t scale_to_decibels(analog_t scale) { return math::log10(scale)*cast<analog_t>(20.0);}
      static gain from_scale(analog_t scale) { return gain(scale_to_decibels(scale)); }
      static gain from_decibels(analog_t db) { return gain(db); }

      // implement casting operators so that gain can be used as a parameter type.
      explicit operator binary_t() const { return db >= 0; }
      explicit operator digital_t() const { return cast<digital_t>(db); }
      explicit operator analog_t() const { return db; }

      analog_t to_scale() const { return decibels_to_scale(db); }
      analog_t to_decibels() const { return db; }
    };
    
    struct duration
    {
      // convert bpm to frequency in Hz
      static constexpr analog_t b_to_f = 1.0 / 60.0;  // NOLINT(clang-diagnostic-implicit-float-conversion)
      static constexpr analog_t f_to_b = 60;

      analog_t samples; // analog so we can express subsample periods.

      duration() : samples(0) {}
      // conversions for parameter
      explicit duration(binary_t b) : samples(b) {}
      explicit duration(size_t s): samples(cast<analog_t>(s)) {}
      explicit duration(digital_t i) : samples(cast<analog_t>(i)) {}
      explicit duration(analog_t a) : samples(a) {}
      explicit operator binary_t() const { return math::abs(samples) >= math::epsilon<analog_t>(); }
      explicit operator digital_t() const { return cast<digital_t>(samples); }
      explicit operator analog_t() const { return samples;}

      static duration from_bpm(analog_t bpm, analog_t sampleRate) { return duration(sampleRate/(bpm*b_to_f)); }
      static duration from_seconds(analog_t seconds, analog_t sampleRate) { return duration(sampleRate*seconds); }
      analog_t to_bpm(analog_t sample_rate) const { return f_to_b*(sample_rate/samples); }
      analog_t to_seconds(analog_t sample_rate) const { return samples/sample_rate; }
    };
  }

  using gain_t     = types::gain;
  using duration_t = types::duration;

  // @todo remove frame namespace, rename channels to frame, move into types namespace
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
      VESSL_INLINE channels() : array<T>(samples, N) { array<T>::fill(0); }
      VESSL_INLINE explicit channels(T m) : array<T>(samples, N) { array<T>::fill(m); }
      VESSL_INLINE channels(const channels& other) : array<T>(samples, N) { other.copyTo(*this); }
      // @todo move and assignment operators

      // mixdown to a mono frame
      channels<T, 1> to_mono() const;
        
      // matrix view of this frame
      math::matrix<T> as_matrix() const;
    };

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator+(channels<T,N> lhs, const channels<T,N>& rhs);

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator-(channels<T,N> lhs, const channels<T,N>& rhs);

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator*(channels<T, N> lhs, const channels<T,N>& rhs);

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator*(channels<T, N> lhs, const T& rhs);

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator*(T lhs, const channels<T, N>& rhs);

    template<typename T, size_t N>
    VESSL_INLINE constexpr channels<T, N> operator^(channels<T, N> lhs, const channels<T,N>& rhs);
  }

  template<typename I, typename O = I>
  class processor
  {
  public:
    using input_t = I;
    using output_t = O;
    
    processor() = default;
    virtual ~processor() = default;
    processor(const processor&) = delete;
    processor(processor&&) = delete;
    processor& operator=(const processor&) = delete;
    processor& operator=(processor&&) = delete;

    virtual output_t process(const input_t& in) = 0;
    virtual void process(source<input_t>& in, sink<output_t>& out);
    virtual void process(array<input_t> in, array<output_t> out);
  };
  
  template<typename T>
  struct process_sample
  {
    process_sample(processor<T>& processor, T sample) : proc_(&processor), sample_(sample) {}
    T& operator>>(T& out) { out = proc_->process(sample_); return out; }
    
  private:
    processor<T>* proc_;
    T sample_;
  };
  
  template<typename T>
  struct process_block
  {
    process_block(processor<T>& processor, array<T> block) : proc_(&processor), block_(block) {}
    // implies in-place processing of block (but what if it didn't? how do?)
    process_block operator>>(processor<T>& next) { proc_->process(block_, block_); return process_block(next, block_); }
    array<T> operator>>(array<T> out) { proc_->process(block_, out); return out; }
    
  private:
    processor<T>* proc_;
    array<T> block_;
  };

  template<typename T>
  process_block<T> operator>>(array<T> in, processor<T>& proc) { return process_block<T>(proc, in); }
  
  template<typename T>
  struct process_stream : source<T>
  {
    process_stream(processor<T>& processor, source<T>& stream) : proc_(&processor), src_(&stream) {}
    binary_t is_empty() const override { return src_->is_empty(); }
    T read() override { return proc_->process(src_->read()); }
    T& operator>>(T& rhs) { rhs = read(); return rhs; }
    process_stream operator>>(processor<T>& processor) { return process_stream(processor, *this); }
    sink<T>& operator>>(sink<T>& out) { proc_->process(*src_, out); return out; }
    array<T> operator>>(array<T> out) { auto w = out.make_writer(); proc_->process(*src_, w); return out; }
    
  private:
    processor<T>* proc_;
    source<T>* src_;
  };
  
  template<typename T>
  process_block<T> operator>>(const T& in, processor<T>& proc) { return process_block<T>(proc, frame::channels<T,1>(in)); }

  template<typename T>
  process_stream<T> operator>>(source<T>& in, processor<T>& proc) { return process_stream<T>(proc, in); }
  
  template<typename T>
  class transform33 : public processor<frame::channels<T,3>>
  {
  public:
    transform33() : matrix_(data_, 3, 3)
    {
      set_identity();
    }

    transform33(const transform33& other) : matrix_(data_, 3, 3)
    {
      memcpy(data_, other.data_, sizeof(T)*3*3);
    }

    [[nodiscard]] VESSL_INLINE math::matrix<T> matrix() const { return matrix_; }
  
    VESSL_INLINE void set_identity() 
    {
      matrix_.clear();
      for (size_t i = 0; i < 3; i++) {
        matrix_.set(i,i, T(1LL));
      }
    }
    
    void set_euler(phase_t pitch, phase_t yaw, phase_t roll);

    VESSL_INLINE void set_euler_radians(analog_t pitchRadians, analog_t yawRadians, analog_t rollRadians)
    {
      return set_euler(cast<phase_t>(pitchRadians / math::two_pi<analog_t>()), 
        cast<phase_t>(yawRadians / math::two_pi<analog_t>()), 
        cast<phase_t>(rollRadians / math::two_pi<analog_t>()));
    }
    
    [[nodiscard]] VESSL_INLINE frame::channels<T,3> process(const frame::channels<T,3>& input) override
    {
      frame::channels<T,3> output;

      output[0] = matrix_[0][0] * input[0] + matrix_[0][1] * input[1] + matrix_[0][2] * input[2];
      output[1] = matrix_[1][0] * input[0] + matrix_[1][1] * input[1] + matrix_[1][2] * input[2];
      output[2] = matrix_[2][0] * input[0] + matrix_[2][1] * input[1] + matrix_[2][2] * input[2];

      // this might be faster?
      //mtrx.multiply(input.toMatrix(), output.toMatrix());

      return output;
    }
    
  private:
    T data_[3 * 3];
    math::matrix<T> matrix_;
  };

  template<typename T>
  class ring : array<T>
  {
  public:
    ring(T* ring_data, size_t data_size) : array<T>(ring_data, data_size), head_(ring_data + data_size - 1) {}

    // expose direct access to underlying array data
    using array<T>::data;
    using array<T>::size;

    void write(const T& v);
    size_t get_write_index() const { return head_ - array<T>::data_; }
    void set_write_index(size_t index) { head_ = array<T>::data_ + index%array<T>::size_; }

    ring operator<<(typename array<T>::reader r);
    
  private:
    T* head_;
  };
  
  class parameter
  {
  public:
    enum class value_type : uint8_t
    {
      none = 0,
      binary = 1, // on/off (binary_t)
      digital = 2, // integral values (digital_t)
      analog = 3, // floating point values (analog_t)
      phase = 4, // phase_t values
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
      value_type    type;
      
      VESSL_INLINE static desc empty() { return {"", 0, value_type::none}; }
    };
    
    template<size_t N>
    struct desc_list
    {
      static constexpr size_t size = N;
      desc list[N];
      
      VESSL_INLINE constexpr desc operator[](id_t id) const
      {
        for (size_t i = 0; i < size; ++i)
        {
          if (list[i].id == id)
          {
            return list[i];
          }
        }
        return desc::empty();
      }
    };
    
    template<typename T>
    struct data
    {
      T value = T(0);
      static constexpr auto type = value_type::user;
    };
    
    template<typename T>
    parameter(const desc& param_desc, const data<T>& param_data) : desc_(param_desc), data_(&const_cast<data<T>&>(param_data).value) {}
      
    // explicitly declared copy-constructors so that they will be used instead of the copy-assignment override
    // when returning parameter objects by value (as in plist::element_at implementations)
    constexpr parameter(parameter& param) : desc_(param.desc_), data_(param.data_) {};
    constexpr parameter(const parameter& param) : desc_(param.desc_), data_(const_cast<void*>(param.data_)) {};
    
    [[nodiscard]] VESSL_INLINE const desc& description() const { return desc_; }
    
    template<typename T>
    VESSL_INLINE T read() const
    {
      switch (desc_.type)
      {
        case value_type::none: return T(0LL);
        case value_type::binary: return cast<T>(*static_cast<binary_t*>(data_));
        case value_type::digital: return cast<T>(*static_cast<digital_t*>(data_));
        case value_type::analog: return cast<T>(*static_cast<analog_t*>(data_));
        case value_type::phase: return cast<T>(*static_cast<phase_t*>(data_));
        case value_type::user: break;
      }
      // attempt to cast data to T, might work?
      return *static_cast<T*>(data_);
    }

    VESSL_INLINE binary_t  read_binary()  const { return read<binary_t>(); }
    VESSL_INLINE digital_t read_digital() const { return read<digital_t>(); }
    VESSL_INLINE analog_t  read_analog()  const { return read<analog_t>(); }
    VESSL_INLINE phase_t   read_phase()   const { return read<phase_t>(); }
    
    VESSL_INLINE explicit operator binary_t()  const { return read<binary_t>(); }
    VESSL_INLINE explicit operator digital_t() const { return read<digital_t>(); }
    VESSL_INLINE explicit operator analog_t()  const { return read<analog_t>(); }
    VESSL_INLINE explicit operator phase_t()   const { return read<phase_t>(); }

    // static-cast T to parameter type before assign
    template<typename T>
    VESSL_INLINE parameter& write(const T& value)
    {
      switch (desc_.type)
      {
        case value_type::none: break;
        case value_type::binary: *static_cast<binary_t*>(data_)   = cast<binary_t>(value); break;
        case value_type::digital: *static_cast<digital_t*>(data_) = cast<digital_t>(value); break;
        case value_type::analog: *static_cast<analog_t*>(data_)   = cast<analog_t>(value); break;
        case value_type::phase: *static_cast<phase_t*>(data_)     = cast<phase_t>(value); break;
        // attempt to cast data to T, might work?
        case value_type::user: *static_cast<T*>(data_) = value; break;
      }
      return *this;
    }
    
    template<typename T>
    VESSL_INLINE parameter& operator=(const T& value)
    {
      write(value);
      return *this;
    }
    
    VESSL_INLINE parameter& operator=(const parameter& rhs)
    {
      if (&rhs != this)
      {
        switch (desc_.type)
        {
        case value_type::none: break;
        case value_type::binary: *static_cast<binary_t*>(data_)   = rhs.read_binary(); break;
        case value_type::digital: *static_cast<digital_t*>(data_) = rhs.read_digital(); break;
        case value_type::analog: *static_cast<analog_t*>(data_)   = rhs.read_analog(); break;
        case value_type::phase: *static_cast<phase_t*>(data_)     = rhs.read_phase(); break;
        case value_type::user: 
          switch (rhs.desc_.type)
          {
          case value_type::none: break;
          case value_type::binary: write(rhs.read_binary()); break;
          case value_type::digital: write(rhs.read_digital()); break;
          case value_type::analog: write(rhs.read_analog()); break;
          case value_type::phase: write(rhs.read_phase()); break;
          case value_type::user: VASSERT(false, "Can't assign user parameter to user parameter with operator="); break;
          }
          break;
        }
      }
      return *this;
    }
    
    static parameter none();
    
  private:
    desc  desc_;
    void* data_;
  };
  
  struct parameter_list : list<parameter>
  {
    // @todo access by ID
  };

  VESSL_INLINE parameter_list::iterator begin(const parameter_list& lst) { return parameter_list::iterator::begin(lst); }
  VESSL_INLINE parameter_list::iterator end(const parameter_list& lst)   { return parameter_list::iterator::end(lst); }
  
  template<size_t N>
  struct plist : parameter_list
  {
    static constexpr size_t num = N;
    VESSL_INLINE size_t size() const override { return N; }
  };
  
  template<>
  struct parameter::data<void*>
  {
    void* value = nullptr;
    static constexpr auto type = value_type::none;
  };
  
  template<>
  struct parameter::data<analog_t>
  {
    analog_t value = 0.f;
    static constexpr auto type = value_type::analog;
  };
  
  template<>
  struct parameter::data<digital_t>
  {
    digital_t value = 0;
    static constexpr auto type = value_type::digital;
  };
  
  template<>
  struct parameter::data<binary_t>
  {
    binary_t value = false;
    static constexpr auto type = value_type::binary;
  };
  
  template<>
  struct parameter::data<phase_t>
  {
    phase_t value = phase_zero;
    static constexpr auto type = value_type::phase;
  };
  
  template<>
  struct parameter::data<gain_t>
  {
    gain_t value = gain_t();
    // @todo add a gain type
    static constexpr auto type = value_type::user;
  };
  
  template<>
  struct parameter::data<duration_t>
  {
    duration_t value = duration_t();
    // @todo add a duration type
    static constexpr auto type = value_type::user;
  };
  
  template<typename T>
  struct param : parameter::data<T>
  {
    parameter operator()(const parameter::desc& d) const { return parameter(d, *this); }
  };
  
  typedef param<analog_t>     analog_p;
  typedef param<digital_t>    digital_p;
  typedef param<binary_t>     binary_p;
  typedef param<phase_t>      phase_p;
  typedef param<gain_t>       gain_p;
  typedef param<duration_t>   duration_p;

  VESSL_INLINE parameter parameter::none() { data<void*> v; return parameter(desc::empty(), v); }
  
  template<typename T>
  process_sample<T> operator>>(const parameter& p, processor<T>& proc) { return process_sample<T>(proc, p.read<T>()); }
  
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
    struct desc
    {
      const char_t*                name;
      const parameter::desc*       params;
      size_t                       param_count;
    };

    virtual ~unit() = default;
    
    // common interface for setting sample rate, for those units that might need it.
    virtual void set_sample_rate(analog_t) {}
    
    // providing a description is optional, but useful.
    virtual desc description() const { return { "", nullptr, 0 }; }
    virtual const parameter_list& parameters() const = 0;
    
    // @todo get parameter by name / id
  };
  
  VESSL_INLINE const parameter::desc* begin(const unit::desc& desc) { return desc.params; }
  VESSL_INLINE const parameter::desc* end(const unit::desc& desc) { return desc.params + desc.param_count; }

  template<typename T>
  class unit_generator : public unit, public generator<T>
  {
  protected:
    explicit unit_generator() : unit(), generator<T>() {}
  };

  template<typename I, typename O = I>
  class unit_processor : public unit, public processor<I,O>
  {
  protected:
    explicit unit_processor() : unit(), processor<I,O>() {}
  };

  namespace interpolation
  {
    template<typename T>
    struct nearest
    {
      T operator()(const T* buffer, analog_t frac_idx);
    };

    template<typename T>
    struct linear
    {
      T operator()(const T* buffer, analog_t frac_idx);
    };

    template<typename T>
    struct cubic
    {
      T operator()(const T* buffer, analog_t frac_idx);
    };
    
    template<typename T, typename I = linear<T>>
    VESSL_INLINE T sample(const T* buffer, analog_t frac_idx)
    {
      VASSERT(frac_idx >= 0, "fracIdx argument to sample must be non-negative");
      static I interpolator;
      return interpolator(buffer, frac_idx);
    }
  };

  namespace random
  {
    // implements the xorshifter algorithm
    static constexpr uint32_t u32_max = UINT32_MAX;
    static uint32_t ru32_seed = 33641;
    VESSL_INLINE void su32(uint32_t seed) {ru32_seed = seed; }
    VESSL_INLINE uint32_t u32()
    {
      ru32_seed ^= ru32_seed << 13; ru32_seed ^= ru32_seed >> 17; ru32_seed ^= ru32_seed << 5;
      return ru32_seed;
    }

    template<typename T>
    VESSL_INLINE T range(T low, T high)
    {
      static constexpr analog_t scale = 1/4294967296.0;
      analog_t r = cast<analog_t>(u32()) * scale;
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
      struct in_out { analog_t operator()(analog_t t) const { return t < 0.5 ? 2*t*t : 1.0 - math::pow<analog_t>(-2*t+2, 2) * 0.5;  } };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
      struct out_in { static out qo; analog_t operator()(analog_t t) const { return t < 0.5 ? qo(2*t) * 0.5f : 1.0 - qo(2*t) * 0.5; } };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
    }

    namespace expo
    {
      struct in { analog_t operator()(analog_t t) const { return t <= math::epsilon<analog_t>() ? 0 : math::pow<analog_t>(2, 10*t-10); } };
      struct out { analog_t operator()(analog_t t) const { return t >= 1.0 - math::epsilon<analog_t>() ? 1.0 : 1.0 - math::pow<analog_t>(2, -10*t);} };  // NOLINT(bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
      struct in_out
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
    VESSL_INLINE T interp(T begin, T end, analog_t t)
    {
      static E ease;
      return (end-begin) * ease(math::constrain(t, 0.f, 1.f)) + begin;
    }

    template<typename T>
    VESSL_INLINE T lerp(T begin, T end, analog_t t) { return interp<linear, T>(begin, end, t); }
    
    // @todo figure out a better name for this
    template<typename T>
    VESSL_INLINE T lerpp(T begin, T end, phase_t t) 
    { 
      return t == phase_zero ? begin 
           : t == phase_360 ? end 
           : begin < end ? begin + (end-begin)*t/phase_360
           : begin - (begin-end)*t/phase_360; 
    }
    
    template<typename T>
    VESSL_INLINE T smooth(T value, T target, analog_t degree = 0.9f) { return value*degree + (1.0 - degree)*target; }
    
    template<>
    VESSL_INLINE digital_t smooth<digital_t>(digital_t value, digital_t target, analog_t degree) { return (value*degree + target)/(degree+1);  }
  }

  // analog unipolar noise generators that generate values in the range [0,1]
  namespace noise
  {
    struct white
    {
      explicit white(analog_t sampleRate) { (void)sampleRate; }
      analog_t operator()() const { return random::range<analog_t>(0, 1); }
    };

    // Implements the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
    // Would be good to dig into the improvements on the algorithm mentioned later in the article.
    struct pink
    {
    public:
      explicit pink(analog_t sample_rate) { (void)sample_rate; }
      
      analog_t operator()()
      {
        int lastKey = key_;
        analog_t sum = 0;
        if (key_ == max_key)
        {
          key_ = 0;
        }
        else
        {
          ++key_;
        }

        int diff = lastKey ^ key_;
        for (int i = 0; i < count; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            white_values_[i] = random::u32() % (range / count);
          }
          sum += cast<analog_t>(white_values_[i]);
        }
        max_sum_ = math::max(sum, max_sum_);
        analog_t n = sum / max_sum_;
        // vassert(!math::isNan(n) && "pink noise generated nan");
        return n;
      }
      
    private:
      static constexpr int range = 128;
      static constexpr int count = 6;
      static constexpr int max_key = 0x1f;

      int key_ = 0;
      analog_t max_sum_ = 90;
      uint32_t white_values_[count] = {
        random::u32() % (range / count), random::u32() % (range / count), random::u32() % (range / count),
        random::u32() % (range / count), random::u32() % (range / count), random::u32() % (range / count)
      };
    };

    // Brownian noise (i.e. random wander) run thru a DC blocking filter.
    // See: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
    // @todo still a bit crunchy
    struct red
    {
      explicit red(analog_t sample_rate) : r_((sample_rate-2.0f)/sample_rate), rc_((1.0f - r_)*200), x_(0), y_(0) {}
      
      analog_t operator()()
      {
        analog_t white = random::range<analog_t>(-rc_, rc_);
        analog_t x = y_ + white;
        // only run the filter when we get close to going out of range
        // to compensate for wandering away from 0.
        y_ = x < -0.49 || x > 0.49f ? x - x_ + r_*y_ : x;
        x_ = x;
        return y_ + 0.5f;
      }

    private:
      analog_t r_;
      analog_t rc_;
      analog_t x_, y_;
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
      gain_t   g;
      args(analog_t sample_rate, analog_t hertz, analog_t kyu, gain_t gain) 
      : sr(sample_rate), hz(hertz), q(kyu), g(gain)
      {
      }
      
      // helper for biquad
      analog_t omega() const { return hz * math::pi<analog_t>() / sr; }
    };
    
    // DC blocking filter, see: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
    template<typename T>
    struct dc_block
    {
      T x1 = T(0), y1 = T(0);
      void process(const T* source, T* dest, size_t block_size, const args& args)
      {
        // gives us about 0.995 for 44100, which is a pretty good R according to the article above.
        analog_t r = (args.sr - 200.0f) / args.sr;
        while (block_size--)
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
      constexpr T butterworth() { return cast<T>(0.70710678118); } // 1/sqrt(2)

      template<typename T>
      constexpr T sallen_key() { return cast<T>(0.5); } 

      template<typename T>
      constexpr T bessel() { return cast<T>(0.57735026919); } // 1/sqrt(3)
    }
    
    template<typename T>
    struct data
    {
      array<analog_t> coeff;
      array<T> state;

      data(analog_t* coeff_data, size_t coeff_size, T* state_data, size_t state_size)
      : coeff(coeff_data, coeff_size), state(state_data, state_size)
      { coeff.fill(0); state.fill(T(0)); }
    };

    // based on https://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
    template<size_t Stages>
    struct biquad
    {
      static constexpr size_t coeff_num = 5;
      
      template<typename T, size_t STATES>
      struct cascade : data<T>
      {
      public:
        cascade() : data<T>(co_, coeff_num * Stages, st_, STATES * Stages), co_{} {}
        
      private:
        analog_t co_[coeff_num*Stages];
        T st_[STATES*Stages];
      };
      
      template<typename T, class CoGen>
      // ReSharper disable once CppInconsistentNaming
      struct df2t final : cascade<T, 2>
      {
        using cascade<T, 2>::coeff;
        using cascade<T, 2>::state;
        using cascade<T, 2>::getCoeffSize;
        using cascade<T, 2>::getStateSize;
        
        static CoGen cg;
         
        void process(const T* source, T* dest, size_t block_size, const args& args);
        // ReSharper disable once CppMemberFunctionMayBeStatic
        size_t stage_count() const { return Stages; }
      };

      template<typename T>
      static void copy(T* coeff);
      
      // coefficient generators
      struct lpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const; };
      struct hpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const; };
      struct bpcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const; };
      struct ntcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const; };
      struct pkcg { void operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t g) const; };
      struct lscg { void operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const; };
      struct hscg { void operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const; };
      
      // base class for filter types to wrap df2T because we specialize it for ARM.
      template<typename T, class CoGen>
      struct flt
      {
        df2t<T, CoGen> df2;
        void process(const T* source, T* dest, size_t block_size, const args& args)
        {
          df2.process(source, dest, block_size, args);
        }
      };

      // filter types for the filter unit generator
      template<typename T>
      struct low_pass final: flt<T, lpcg> {};
      
      template<typename T>
      struct high_pass final : flt<T, hpcg> {};
      
      template<typename T>
      struct band_pass final : flt<T, bpcg> {};
        
      template<typename T>
      struct notch final : flt<T, ntcg> {};
        
      template<typename T>
      struct peak final : flt<T, pkcg> {};
        
      template<typename T>
      struct low_shelf final : flt<T, lscg> {};
        
      template<typename T>
      struct high_shelf final : flt<T, hscg> {};
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
  
  // classes can subclass this to add support for tempo detection of a clock signal (i.e. pulse train)
  class clockable
  {
  public:
    using period_t = uint32_t;

    clockable(analog_t sample_rate, period_t sample_period_min, period_t sample_period_max, analog_t bpm = 60)
    : tempo_(duration_t::from_bpm(bpm, sample_rate)), period_min_(sample_period_min), period_max_(sample_period_max)
    , ticks_(0), sample_rate_(sample_rate) {}
    virtual ~clockable() = default;
    clockable(const clockable&) = default;
    clockable(clockable&&) = default;
    clockable& operator=(const clockable&) = default;
    clockable& operator=(clockable&&) = default;
    
    // users should call clock at the beginning of every clock pulse
    void clock();
    void clock(period_t sample_delay);

    analog_t bpm() const { return tempo_.to_bpm(sample_rate_); }
    // length of one clock pulse in samples
    analog_t period() const { return tempo_.samples;}
    
  protected:
    // subclasses should call tick for every sample generated/processed
    void tick() { ++ticks_; }
    void tick(size_t t) { ticks_ += t; }

    // subclasses can override this to be notified every time they receive a clock pulse
    virtual void tock(size_t sample_delay) { (void)sample_delay; }
    
    duration_t tempo_;
    period_t period_min_;
    period_t period_max_;
    period_t ticks_;
    analog_t sample_rate_;
  };
  
  template<typename T>
  struct smoother
  {
    T        value;
    analog_t degree;
    
    explicit smoother(analog_t smoothing_degree = 0.9f, T initial_value = T(0))
      : value(initial_value), degree(smoothing_degree)
    {
    }
    
    VESSL_INLINE explicit operator T() const { return value; }
    
    // so we can use this like OWL's SmoothValue
    VESSL_INLINE T operator=(const T& v)
    {
      analog_t d = math::constrain<analog_t>(degree, 0.0, 1.0);
      return value = easing::smooth(value, v, d);
    }
    
    VESSL_INLINE T operator=(const parameter& p)
    {
      return *this = p.read<T>();
    }
  };
  
  // a waveform that can be evaluated using a normalized phase value
  // implementors should accept negative phase, as well as phase values outside [-1,1]
  template<typename T>
  struct waveform
  {
    using sample_t = T;
    
    waveform() = default;
    virtual ~waveform() = default;
    waveform(const waveform&) = default;
    waveform(waveform&&) = default;
    waveform& operator=(const waveform&) = default;
    waveform& operator=(waveform&&) = default;

    virtual sample_t evaluate(phase_t phase) const = 0;  // NOLINT(portability-template-virtual-member-function)
  };

  namespace waves
  {
    template<typename T>
    struct sine final : waveform<T>
    {
      // default implementation assumes default parameter type (floating point)
      VESSL_INLINE T evaluate(phase_t phase) const override { return math::sin<T>(phase); }
    };

    template<typename T>
    struct cosine final : waveform<T>
    {
      VESSL_INLINE T evaluate(phase_t phase) const override { return math::cos<T>(phase); }
    };

    template<typename T>
    struct square final : waveform<T>
    {
      phase_t pulse_width;
      square() : pulse_width(phase_180) {}
      explicit square(phase_t pw) : pulse_width(pw) {}
      VESSL_INLINE T evaluate(phase_t phase) const override { return phase < pulse_width ? 1 : -1; }
    };

    // same as square, but unipolar
    template<typename T>
    struct clock final : waveform<T>
    {
      phase_t pulse_width;
      clock() : pulse_width(phase_180) {}
      explicit clock(phase_t pw) : pulse_width(pw) {}
      VESSL_INLINE T evaluate(phase_t phase) const override { return phase < pulse_width ? 1 : 0; }
    };
  }

  // a fixed-sized buffer that supports sampling it a normalized phase.
  // both positive and negative phases are supported
  template<typename T, size_t N, typename I = interpolation::linear<T>>
  class wavetable final : public waveform<T>
  {
  public:
    wavetable(): waveform<T>() {}

    // assumes source can provide at least N samples
    explicit wavetable(source<T>& source);
    explicit wavetable(const waveform<T>& waveform);

    // ReSharper disable once CppMemberFunctionMayBeStatic
    VESSL_INLINE size_t size() const { return N; }

    VESSL_INLINE T get(const size_t i) const { return buffer[i + 1]; }
    VESSL_INLINE void set(const size_t i, T val);
    
    // implement waveform:
    T evaluate(phase_t phase) const override;
    
  private:
    T buffer[N + 3] = {};
  };
  
  // generates stepped analog noise at rate values per second.
  template<typename T, typename N>
  class noise_generator : public unit_generator<T>, protected plist<1>
  {
  public:
    N noise_source;
    
    explicit noise_generator(analog_t sample_rate = 1)
    : unit_generator<T>(), noise_source(sample_rate), dt_(1.0f/sample_rate), step_(0)
    {
      value_ = noise_source();
      next_ = noise_source();
      params_.rate.value = sample_rate;
    }
    
    void set_sample_rate(float sample_rate) override { dt_ = 1.0f / sample_rate;}
    const parameter_list& parameters() const override { return *this; }

    parameter rate() const { return params_.rate({ "rate", 'r', analog_p::type }); }

    // generates stepped noise in the range [0,1] at the given rate
    T generate() override;

    // smooths the stepped noise with the given easing
    template<typename E>
    T generate();

  protected:
    parameter element_at(size_t index) const override { parameter p[num] = { rate() }; return p[index]; }

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
    parameter from() const { return params_.from({ "from", 'f', param<T>::type }); }
    parameter to() const { return params_.to({ "to", 't', param<T>::type }); }
    parameter duration() const { return params_.duration({ "duration", 'd', analog_p::type }); }

    // outs
    parameter eor() const { return params_.eor({ "eor", 'e', binary_p::type }); }
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

      parameter target() const { return params_.target({ "target", 't', analog_p::type }); }
      parameter duration() const { return params_.duration({ "duration", 'd', analog_p::type }); }
      parameter active() const { return params_.active({ "active", 'a', binary_p::type }); }
      parameter eos() const { return params_.eos({ "eos", 'e', binary_p::type }); }
      // current value of the stage
      parameter value() const { return params_.output({ "value", 'v', analog_p::type }); }

      void start(T from_value) { begin_ = from_value; params_.output.value = from_value; time_ = -dt_; params_.active.value = true; params_.eos.value = false; }
      void reset() { params_.active.value = false; params_.eos.value = false; time_ = -dt_; params_.output.value = 0; }
      
      template<typename E>
      T generate() { return active() ? step<E>() : params_.output.value; }
      T generate() override { return generate<easing::linear>(); }
      
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
    
    void set_sample_rate(float sampleRate) override;
    
    stage& get_stage(size_t idx) { return idx == stages_.size() ? final_ : stages_[idx]; }
    const stage& get_stage(size_t idx) const { return idx == stages_.size() ? final_ : stages_[idx]; }
    size_t get_stage_count() const { return stages_.size() + 1; }

    stage& current_stage() { return get_stage(stage_idx_); }
    const stage& current_stage() const { return get_stage(stage_idx_);}
    stage& final_stage() { return final_; }
    const stage& final_stage() const { return final_; }
    
    const parameter_list& parameters() const override { return *this; }

    parameter value() const { return current_stage().value(); }
    parameter eoc() const { return params_.eoc({ "eoc", 'e', binary_p::type }); }

    // make this a parameter we check in generate?
    virtual void trigger();

    template<typename E>
    T generate();
    T generate() override { return generate<easing::linear>(); }

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
    typename envelope<T>::stage& decay() { return envelope<T>::finalStage(); }
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
    typename envelope<T>::stage& release() { return envelope<T>::finalStage(); }
    using ad<T>::eoc; 

    void gate(T value);
    void gate(binary_t on) { gate(on ? T(1) : T(0)); }
    void trigger() override { attack().target() = T(1); ad<T>::trigger(); gate_on_ = false; }
    using envelope<T>::generate;

  protected:
    binary_t should_advance(size_t current_stage_idx) override
    {
      return ad<T>::shouldAdvance(current_stage_idx) && !gate_on_;
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
      envelope<T>::finalStage().duration() = release_duration;
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
      return current_stage_idx == 1 ? decay_stage_.eos() && !gate_on_ : envelope<T>::shouldAdvance(current_stage_idx);
    }
    
  private:
    typename envelope<T>::stage attack_stage_;
    typename envelope<T>::stage decay_stage_;
    binary_t gate_on_;
  };
  
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

    parameter rise() const { return params_.rise({ "rise", 'a', analog_p::type }); }
    parameter fall() const { return params_.fall({ "fall", 'd', analog_p::type }); }
    parameter rising() const { return params_.rising({ "rising", 'r', binary_p::type }); }
    parameter falling() const { return params_.falling({ "falling", 'f', binary_p::type }); }
    parameter value() const { return params_.output({ "value", 'v', param<T>::type }); }
    
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

  // note: W must implement waveform<T>
  template<class W>
  class oscil final : public unit_generator<typename W::sample_t>, protected plist<4>
  {
  public:
    W waveform;
    using sample_t = typename W::sample_t;
    
    oscil() : unit_generator<sample_t>(), phase_(phase_zero), dt_(phase_zero) { params_.fHz.value = 440.0; }

    template<typename... Ts>
    explicit oscil(analog_t sample_rate, analog_t freq_in_hz, Ts... wargs)
    : unit_generator<sample_t>(), waveform(wargs...), phase_(phase_zero), dt_(cast<phase_t>(1.0f/sample_rate))
    { params_.fHz.value = freq_in_hz; }
    
    void set_sample_rate(float sample_rate) override { dt_ = cast<phase_t>(1.0f/sample_rate); }
    
    [[nodiscard]] const parameter_list& parameters() const override { return *this; }
    
    // frequency in Hz without FM applied
    [[nodiscard]] parameter fhz() const { return params_.fHz({ "frequency", 'f', analog_p::type }); }
    // linear frequency modulation
    [[nodiscard]] parameter fm_lin() const { return params_.fmLin({ "fm (lin)", 'l', analog_p::type }); }
    // v/oct (exponential) frequency modulation
    [[nodiscard]] parameter fm_exp() const { return params_.fmExp({ "fm (v/oct)", 'v', analog_p::type }); }
    // phase modulation
    [[nodiscard]] parameter pm() const { return params_.pm({ "phase mod", 'p', phase_p::type }); }

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
  class delay_line final : public ring<T>, public waveform<T>
  {
  public:
    delay_line(T* delay_line_data, size_t data_size) : ring<T>(delay_line_data, data_size) {}

    using ring<T>::data;
    using ring<T>::size;
    using ring<T>::get_write_index;
    using ring<T>::set_write_index;

    // reads behind the write head with sampleDelay (i.e. the ith sample previously written)
    // where a delay of 0 samples will give the most recently written value.
    T read(size_t sample_delay) const;

    // reads behind the write head with a fractional sampleDelay and given interpolation
    template<typename I>
    T read(analog_t sample_delay) const;

    // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
    T evaluate(phase_t phase) const override;
  };
  
  template<typename T, typename I = interpolation::linear<T>>
  class delay : public unit_processor<T>, protected plist<2>
  {
  public:
    delay(array<T> delay_buffer, analog_t sample_rate, analog_t delay_in_seconds = 0, analog_t feedback_amount = 0)
    : unit_processor<T>(), buffer_(delay_buffer.getData(), delay_buffer.size()), dt_(1.0f/sample_rate)
    {
      params_.time.value = duration_t::from_seconds(delay_in_seconds, sample_rate);
      delay_in_samples_ = params_.time.value.samples;
      params_.feedback.value = feedback_amount;
    }
    
    void set_sample_rate(float sampleRate) override { dt_ = 1.0f / sampleRate; }
    const parameter_list& parameters() const override { return *this; }

    delay_line<T>& buffer() { return buffer_; }
    const delay_line<T>& buffer() const { return buffer_; }

    /// delay time expressed as vessl::duration (i.e. samples), can be set using an analog_t
    parameter time() const { return params_.time({ "time", 't', duration_p::type }); }
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    parameter feedback() const { return params_.feedback({ "feedback", 'f', analog_p::type }); }
    
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
    
    delay_line<T> buffer_;
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
    
    parameter response() const { return params_.response({ "response time", 'r', analog_p::type }); }

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
    explicit freeze(array<T> freeze_buffer, float sample_rate) : unit()
    , delay_line_(freeze_buffer.getData(), freeze_buffer.size())
    , phase_(0), crossfade_(0.75f)
    , freeze_delay_(0), freeze_size_(0)
    , read_rate_(1), dt_(1.0/sample_rate)
    {
      params_.size.value.samples = delay_line_.size()-1;
      freeze_size_ = params_.size.value.samples;
      params_.rate.value = 1.0;
    }
    
    void set_sample_rate(analog_t sr) override { dt_ = 1.0f/sr; }
    const parameter_list& parameters() const override { return *this; }
    
    delay_line<T>& get_delay_line() { return delay_line_; }
    const delay_line<T>& get_delay_line() const { return delay_line_; }
    
    parameter enabled() const { return params_.enabled({ "enabled", 'e', binary_p::type }); }
    // end of the freeze loop in samples relative to the most recently recorded sample
    parameter position() const { return params_.position({ "position", 'p', analog_p::type }); }
    // size of the freeze loop as a duration (samples).
    // the beginning of the freeze loop, when played forward, will be position + size.
    parameter duration() const { return params_.size({ "duration", 'd', analog_p::type }); }
    // rate of playback when enabled, can be negative to play in reverse
    parameter rate() const { return params_.rate({ "rate", 'r', analog_p::type }); }
    // should this be a parameter? there's not much gained by it.
    analog_t phase() const { return phase_; }
    // reset the phase to zero (argument for making it a parameter?)
    void reset() { phase_ = 0; }

    T generate() override;

    T process(const T& in) override;

    void process(source<T>& source, sink<T>& sink) override { processor<T>::process(source, sink); }
    
    template<time::mode TimeMode = time::mode::slew>
    void generate(array<T> output) { proc_gen<TimeMode, false>(output, output); }

    template<time::mode TimeMode = time::mode::slew>
    void process(array<T> input, array<T> output) { proc_gen<TimeMode, true>(input, output); }
    
  protected:
    parameter element_at(size_t index) const override
    {
      parameter p[num] = { enabled(), position(), size(), rate() };
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
    
    delay_line<T> delay_line_;
    analog_t phase_;
    // used to crossfade between the incoming signal and the freeze signal when enabled changes.
    smoother<analog_t> crossfade_;
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
      params_.fHz.value = freq_in_hz;
      params_.q.value = kyu;
      params_.emphasis.value = emphasis;
    }
    
    void set_sample_rate(analog_t sample_rate) override { sample_rate_ = sample_rate; }
    const parameter_list& parameters() const override { return *this; }

    parameter fhz() const { return params_.fhz({ "fHz", 'f', analog_p::type }); }
    parameter q() const { return params_.q({ "q", 'q', analog_p::type }); }
    // unused by some filter types (see filtering section)
    parameter emphasis() const { return params_.emphasis({ "emphasis", 'e', analog_p::type }); }
    
    T process(const T& in) override
    {
      T out;
      filtering::args fargs = {
        sample_rate_,
        params_.fHz.value,
        math::max(params_.q.value, 0.01),
        params_.emphasis.value
      };
      func_.process(&in, &out, 1, fargs);
      return out;
    }
    
    void process(array<T> in, array<T> out) override
    {
      filtering::args fargs = {
        sample_rate_,
        params_.fHz.value,
        math::max(params_.q.value, 0.01),
        params_.emphasis.value
      };
      func_.process(in.getData(), out.getData(), in.size(), fargs);
    }
    
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

    parameter rate() const { return params_.bitRate({ "bit rate", 'r', analog_p::type }); }
    parameter depth() const { return params_.bitDepth({ "bit depth", 'd', analog_p::type }); }
    parameter mangle() const { return params_.mangle({ "mangle", 'm', binary_p::type }); }

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

  // A simple peak limiter adapted from pinchenettes/stmlib via DaisySP
  template<typename T>
  class limiter : public unit_processor<T>, protected plist<2>
  {
  public:
    explicit limiter(gain_t pre_gain = gain_t::from_decibels(0)) : unit_processor<T>()
    {
      params_.preGain.value = pre_gain;
      params_.peak.value = 0.5;
    }
    
    const parameter_list& parameters() const override { return *this; }

    parameter pre_gain() const { return params_.preGain({ "pre-gain", 'g', gain_p::type }); }
    parameter peak() const { return params_.peak({ "peak", 'k', param<T>::type }); }

    T process(const T& in) override;
    using unit_processor<T>::process;
    
  protected:
    parameter element_at(size_t index) const override
    {
      parameter p[num] = { pre_gain(), peak() };
      return p[index];
    }
    
  private:
    struct
    {
      gain_p pre_gain;
      param<T> peak;
    } params_;
  };
}
  
/////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation
//

#ifdef ARM_CORTEX
#include "vessl_arm_math.inl"
#endif

#include "vessl_qmath.inl"
#include "vessl_frame.inl"

namespace vessl
{
  namespace  math
  {
    template<>
    VESSL_INLINE analog_t sin<analog_t, phase_t>(phase_t z) 
    { 
      return math::sin<analog_t>(math::two_pi<analog_t>() * cast<analog_t>(z)); 
    }

    template<>
    VESSL_INLINE analog_t cos<analog_t, phase_t>(phase_t z) 
    { 
      return math::cos<analog_t>(math::two_pi<analog_t>() * cast<analog_t>(z)); 
    }
  
    template<typename T, size_t N>
    VESSL_INLINE frame::channels<T, N> round(frame::channels<T, N> x)
    {
      frame::channels<T, N> result;
      for (size_t i = 0; i < N; i++)
      {
        result[i] = round(x[i]);
      }
      return result;
    }
  }
  
  template<typename I, typename O>
  void processor<I,O>::process(source<I>& in, sink<O>& out)
  {
    while (in && out)
    {
      out.write(process(in.read()));
    }
  }

  template<typename I, typename O>
  void processor<I,O>::process(array<I> in, array<O> out)
  {
    auto r = in.make_reader();
    auto w = out.make_writer();
    process(r, w);
  }

  // @todo ARM specialization
  template<typename T>
  void array<T>::writer::write(const reader& r)
  {
    size_t rsz = r.available();
    VASSERT(available() >= rsz, "Not enough space in writer for the contents of reader");
    const T* rh = *r;
    memcpy(static_cast<void*>(head_), static_cast<const void*>(rh), rsz * sizeof(T));
    head_ += rsz;
  }
  
  template<typename T>
  void array<T>::copy_to(array dest)
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
    VASSERT(size_ <= dest.size, "arrays are have different lengths or destination is too small");
    reader a(data_, size_);
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
    VASSERT(size_ == other.size_
      && size_ <= dest.size_
      , "arrays are have different lengths or destination is too small");
    reader a(data_, size_);
    reader b(other.data_, other.size_);
    writer c(dest.data_, dest.size_);
    while (a)
    {
      c << a.read() + b.read();
    }
    return dest;
  }

  template<typename T>
  array<T> array<T>::subtract(array other, array dest) const
  {
    VASSERT(size_ == other.size_ 
      && size_ <= dest.size_
      , "arrays are have different lengths or destination is too small");
    reader a(data_, size_);
    reader b(other.data_, other.size_);
    writer c(dest.data_, dest.size_);
    while (a)
    {
      c << a.read() - b.read();
    }
    return dest;
  }

  template<typename T>
  array<T> array<T>::scale(T value, array dest) const
  {
    VASSERT(size_ <= dest.size_, "destination size is too small");
    reader a(data_, size_);
    writer b(dest.data_, dest.size_);
    while (a)
    {
      b << a.read() * value;
    }
    return dest;
  }

  template<typename T>
  array<T> array<T>::multiply(array other, array dest) const
  {
    VASSERT(size_ == other.size_
      && size_ <= dest.size_
      , "arrays are have different lengths or destination is too small");
    reader a(data_, size_);
    reader b(other.data_, other.size_);
    writer c(dest.data_, dest.size_);
    while (a)
    {
      c << a.read() * b.read();
    }
    return dest;
  }

  template<typename T>
  typename array<T>::writer& operator<<(typename array<T>::writer& w, const T& v)
  {
    w.write(v);
    return w;
  }

  namespace math
  {
  template<typename T>
  matrix<T> matrix<T>::add(matrix other, matrix dest) const
  {
    VASSERT(rows() == other.rows() 
      && rows() == dest.rows() 
      && columns() == other.columns() 
      && columns() == dest.columns()
      , "matrices do not have the same dimensions");
    
    array<T> lhs(data(), size());
    array<T> rhs(other.data(), other.size());
    array<T> dst(dest.data(), dest.size());
    lhs.add(rhs, dst);
    return dest;
  }

  template<typename T>
  matrix<T> matrix<T>::subtract(matrix other, matrix dest) const
  {
    VASSERT(rows() == other.rows() 
      && rows() == dest.rows() 
      && columns() == other.columns() 
      && columns() == dest.columns()
      , "matrices do not have the same dimensions");
    
    array<T> lhs(data(), size());
    array<T> rhs(other.data(), other.size());
    array<T> dst(dest.data(), dest.size());
    lhs.subtract(rhs, dst);
    return dest;
  }
  
  template<typename T>
  matrix<T> matrix<T>::scale(T value, matrix dest) const
  {
    array<T> lhs(data(), size());
    array<T> dst(dest.data(), dest.size());
    lhs.scale(value, dst);
    return dest;
  }
  
  template<typename T>
  matrix<T> matrix<T>::multiply(matrix other, matrix dest) const
  {
    VASSERT(columns() == other.rows(), "Incompatible matrix sizes in operands");
    VASSERT(dest.rows() == rows(), "Incorrect number of rows in destination");
    VASSERT(dest.columns() == other.columns(), "Incorrect number of columns in destination");
    
    for(size_t i = 0; i < rows(); i++)
    {
      for(size_t j = 0; j < other.columns(); j++)
      {
        T accum = T(0LL);
        for(size_t k = 0; k < other.rows(); k++)
        {
          accum += get(i, k) * other.get(k, j);
        }
        dest.set(i, j, accum);
      }
    }
    return dest;
  }
  
  template<typename T>
  array<T> matrix<T>::multiply(const array<T>& vector, array<T> dest) const
  {
    VASSERT(columns() == vector.size(), "Incompatible operands");
    VASSERT(dest.size() == rows(), "Incompatible destination size");
    
    for(size_t i = 0; i < rows(); i++)
    {
      T accum = 0;
      for(size_t j = 0; j < vector.size(); j++)
      {
        accum += get(i, j) * vector[j];
      }
      dest[i] = accum;
    }
    return dest;
  }
  }

  template <typename T>
  void transform33<T>::set_euler(phase_t pitch, phase_t yaw, phase_t roll)
  {
    T cosa = math::cos<T>(roll);
    T sina = math::sin<T>(roll);

    T cosb = math::cos<T>(yaw);
    T sinb = math::sin<T>(yaw);

    T cosc = math::cos<T>(pitch);
    T sinc = math::sin<T>(pitch);

    // row*3 + col
    matrix_.set(0, 0, cosa * cosb);
    matrix_.set(0, 1, cosa * sinb*sinc - sina*cosc);
    matrix_.set(0, 2, cosa * sinb*cosc + sina*sinc);

    matrix_.set(1, 0, sina * cosb);
    matrix_.set(1, 1, sina * sinb*sinc + cosa*cosc);
    matrix_.set(1, 2, sina * sinb*cosc - cosa*sinc);

    matrix_.set(2, 0, -sinb);
    matrix_.set(2, 1, cosb * sinc);
    matrix_.set(2, 2, cosb * cosc);
  }

  template<typename T>
  void ring<T>::write(const T& v)
  {
    *head_++ = v;
    if (head_ == array<T>::end())
    {
      head_ = array<T>::begin();
    }
  }

  template<typename T>
  ring<T> ring<T>::operator<<(typename array<T>::reader r)
  {
    VASSERT(r.available() < size(), "reader size is larger than ring size");
    while (r)
    {
      write(r.read());
    }
    return *this;
  }

  inline void clockable::clock()
  {
    tempo_.samples = cast<analog_t>(math::constrain(ticks_, period_min_, period_max_));
    ticks_ = 0;
    tock(0);
  }

  inline void clockable::clock(period_t sample_delay)
  {
    tempo_.samples = cast<analog_t>(math::constrain(ticks_ + sample_delay, period_min_, period_max_));
    ticks_ = 0;
    tock(sample_delay);
  }

  template<typename T>
  T interpolation::nearest<T>::operator()(const T* buffer, analog_t frac_idx)
  {
    return buffer[cast<size_t>(math::round(frac_idx))];
  }

  template<typename T>
  T interpolation::linear<T>::operator()(const T* buffer, analog_t frac_idx)
  {
    analog_t idx;
    analog_t frac = math::mod(frac_idx, &idx);
    size_t x0 = cast<size_t>(idx);
    return buffer[x0] + (buffer[x0 + 1] - buffer[x0]) * frac;
  }

  template<typename T>
  T interpolation::cubic<T>::operator()(const T* buffer, analog_t frac_idx)
  {
    static constexpr analog_t div6 = (1. / 6.);
    static constexpr analog_t div2 = (0.5);

    analog_t idx;
    analog_t frc = math::mod(frac_idx, &idx);
    analog_t fm1 = frc - 1.f;
    analog_t fm2 = frc - 2.f;
    analog_t fp1 = frc + 1.f;
    size_t x0 = idx;
    return -frc * fm1 * fm2 * div6 * buffer[x0 - 1] 
          + fp1 * fm1 * fm2 * div2 * buffer[x0] 
          - fp1 * frc * fm2 * div2 * buffer[x0 + 1] 
          + fp1 * frc * fm1 * div6 * buffer[x0 + 2];
  }
  
  template<size_t Stages>
  template<typename T, class CoGen>
  void filtering::biquad<Stages>::df2t<T, CoGen>::process(const T* source, T* dest, size_t block_size, const args& args)
  {
    // update our coefficients
    cg(coeff.data(), args.omega(), args.q, args.g);
    
    // run the filter
    const T* input = source;
    array<analog_t>::reader cor = coeff.getReader();
    for (size_t s = 0; s < Stages; s++)
    {
      analog_t b0 = cor.read(); analog_t b1 = cor.read(); analog_t b2 = cor.read();
      analog_t a1 = cor.read();  analog_t a2 = cor.read();
      T* st = state.data() + 2*s;
      T d1 = st[0]; T d2 = st[1];
      typename array<T>::reader r(input, block_size);
      typename array<T>::writer w(dest, block_size);
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

  template<size_t Stages>
  template<typename T>
  void filtering::biquad<Stages>::copy(T* coeff) 
  {
    if (Stages > 1)
    {
      array<T> src(coeff, coeff_num);
      for (size_t i = 1; i < Stages; i++)
      {
        array<T> dst(coeff + coeff_num*i, coeff_num);
        src.copy_to(dst);
      }
    }
  }

  template<size_t Stages>
  void filtering::biquad<Stages>::lpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
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

  template<size_t Stages>
  void filtering::biquad<Stages>::hpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
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

  template<size_t Stages>
  void filtering::biquad<Stages>::bpcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
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

  template<size_t Stages>
  void filtering::biquad<Stages>::ntcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t _) const
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

  template<size_t Stages>
  void filtering::biquad<Stages>::pkcg::operator()(analog_t* coeff, analog_t omega, analog_t q, gain_t g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.to_decibels())/20);
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

  template<size_t Stages>
  void filtering::biquad<Stages>::lscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.to_decibels())/20);
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

  template<size_t Stages>
  void filtering::biquad<Stages>::hscg::operator()(analog_t* coeff, analog_t omega, analog_t _, gain_t g) const
  {
    analog_t K = math::tan(omega);
    analog_t V = math::exp10(math::abs(g.to_decibels())/20);
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
  T wavetable<T, N, I>::evaluate(phase_t phase) const
  {
    analog_t idx = cast<analog_t>(phase) * N;
    return interpolation::sample<T, I>(buffer, idx);
  }

  template<typename T, typename N>
  T noise_generator<T, N>::generate()
  {
    step_ += dt_ * rate();
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
  T noise_generator<T, N>::generate()
  {
    T s = generate();
    return easing::interp<E>(s, next_, step_);
  }

  template<typename T>
  void ramp<T>::trigger()
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
  T ramp<T>::generate()
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
    return easing::lerp(params_.from.value, params_.to.value, lt);
  }

  template<typename T>
  template<typename E>
  T envelope<T>::stage::step()
  {
    analog_t s = dt_ / math::max<analog_t>(cast<analog_t>(duration()), dt_);
    time_ += s;
    analog_t t = math::constrain<analog_t>(time_, 0.0, 1.0);
    params_.output.value = easing::interp<E, T>(begin_, params_.target.value, t);
    if (time_ >= 1)
    {
      params_.active.value = false;
      params_.eos.value = true;
    }
    return params_.output.value;
  }

  template<typename T>
  void envelope<T>::trigger()
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
  T envelope<T>::generate()
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
  void envelope<T>::set_sample_rate(analog_t sample_rate)
  {
    for (stage& s : stages_)
    {
      s.set_sample_rate(sample_rate);
    }
    final_.setSampleRate(sample_rate);
  }

  template<typename T>
  void asr<T>::gate(T value)
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

  template<typename T>
  T slew<T>::process(const T& v)
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

  template<class W>
  VESSL_INLINE typename W::sample_t oscil<W>::generate()
  {
    typename W::sample_t val = waveform.evaluate(phase_ + params_.pm.value);
    analog_t f = params_.fHz.value * math::exp2(params_.fmExp.value) + params_.fmLin.value;
    phase_ += static_cast<phase_t>(dt_ * f);
    return val;
  }

  template <class W>
  VESSL_INLINE void oscil<W>::generate(sink<typename W::sample_t>& dest)
  {
    analog_t freq = params_.fHz.value * math::exp2(params_.fmExp.value) + params_.fmLin.value;
    phase_t step = static_cast<phase_t>(dt_ * freq);
    while(!dest.isFull())
    {
      dest << waveform.evaluate(phase_ + params_.pm.value);
      phase_ += step;
    }
  }

  template<typename T>
  T delay_line<T>::read(size_t sample_delay) const
  {
    assert(sample_delay < size());
    sample_delay = size() - 1 - sample_delay;
    size_t idx = get_write_index() + sample_delay;
    return data()[idx % size()];
  }

  template<typename T>
  template<typename I>
  T delay_line<T>::read(analog_t sample_delay) const
  {
    assert(sample_delay >= 0 && sample_delay < size());
    analog_t fSize = cast<analog_t>(size());
    sample_delay = fSize - 1 - sample_delay;
    analog_t fidx = cast<analog_t>(get_write_index()) + sample_delay;
    analog_t idx;
    analog_t f = math::mod(fidx, &idx);
    size_t x0 = cast<size_t>(idx) % size();
    size_t x1 = (x0 + 1) % size();
    size_t x2 = (x0 + 2) % size();
    const T* d = data();
    T s[3] = { d[x0], d[x1], d[x2] };
    return interpolation::sample<T, I>(s, f);
  }

  template<typename T>
  T delay_line<T>::evaluate(phase_t phase) const
  {
    analog_t fSize = cast<analog_t>(size());
    analog_t sampleDelay = cast<analog_t>(phase_360 - phase) * fSize;
    return read<interpolation::linear<T>>(sampleDelay);
  }

  template<typename T, typename I>
  T delay<T, I>::process(const T& in)
  {
    delay_in_samples_ = params_.time.value.samples;
    // delay time in samples
    analog_t dts = math::constrain<analog_t>(delay_in_samples_, 0.f, cast<analog_t>(buffer_.size()-1));
    analog_t s = buffer_.template read<I>(dts);
    analog_t fbk = math::constrain<analog_t>(cast<analog_t>(feedback()), -1.0, 1.0);
    buffer_.write(in + s * fbk);
    return s;
  }

  template<typename T, typename I>
  template<time::mode TimeMode>
  void delay<T, I>::process(array<T> input, array<T> output)
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
        delay_in_samples_ = easing::lerp(delay_in_samples_, params_.time.value.samples, dst);
        // delay time in samples
        analog_t dts = math::constrain<analog_t>(delay_in_samples_, 0.f, cast<analog_t>(buffer_.size()-1));
        analog_t wet = buffer_.template read<I>(dts);
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
        T wet = (1.0f - fade) * buffer_.template read<I>(fts) + fade * buffer_.template read<I>(tts);
        buffer_.write(in + wet*fbk);
        w << wet;
        fade += fadeInc;
      }
        
      delay_in_samples_ = targetSampleDelay;
    }
  }

  template<typename T>
  T follow<T>::process(const T& in) 
  {
    writer_.write(in);
    if (writer_.isFull())
    {
      previous_ = current_;
      current_ = T(0);
      auto r = window_.make_reader();
      while (r)
      {
        current_ *= delta_;
        current_ += (1.0 - delta_)*math::abs(r.read());
      }
      writer_ = window_.writer();
    }

    analog_t t = 1.0 - cast<analog_t>(writer_.available()) / cast<analog_t>(window_.size());
    return previous_ + (current_ - previous_) * t;
  }

  template<typename T, typename I>
  T freeze<T, I>::generate() 
  {
    freeze_delay_ = cast<analog_t>(position());
    freeze_size_  = params_.size.value.samples;
    analog_t sampleDelay = freeze_delay_ + (1.0-phase_)*freeze_size_;
    phase_ = math::wrap01(phase_ + rate() / freeze_size_);
    return delay_line_.template read<I>(sampleDelay);
  }

  template<typename T, typename I>
  T freeze<T, I>::process(const T& in) 
  {
    binary_t is_enabled = cast<binary_t>(enabled());
    analog_t wet_level = crossfade_ = (is_enabled ? 1.0 : 0.0);
    T wet = generate();
    if (!is_enabled)
    {
      delay_line_.write(in);
    }
    return wet_level*wet + (1.0 - wet_level)*in;
  }

  template<typename T, typename I>
  template<time::mode TimeMode, bool UseInput>
  void freeze<T, I>::proc_gen(array<T> input, array<T> output) 
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
      analog_t fs = params_.size.value.samples;
      analog_t rt = params_.rate.value;
      analog_t st = dt_;
      while (w)
      {
        freeze_delay_ = easing::lerp(freeze_delay_, fd, st*20);
        freeze_size_ = easing::lerp(freeze_size_, fs, st*20);
        read_rate_ = easing::lerp(read_rate_, rt, st*10);
        analog_t sample_delay = freeze_delay_ + (1.0 - phase_)*freeze_size_;
        phase_ = math::wrap01(phase_ + read_rate_/freeze_size_);
        T wet = delay_line_.template read<I>(sample_delay);

        if (UseInput)
        {
          binary_t is_enabled = params_.enabled.value;
          analog_t wet_level = crossfade_ = (is_enabled ? 1.0 : 0.0);
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
      analog_t fs0 = freeze_size_,    fs1 = params_.size.value.samples;
      analog_t fade = 0, fadeInc = 1.0 / output.size();
      analog_t r0 = read_rate_, r1 = params_.rate.value;
      analog_t p0 = phase_, dp0 = r0/fs0, dp1 = r1/fs1;
      while (w)
      {
        analog_t sd0 = fd0 + fs0*(1.0-p0);
        analog_t sd1 = fd1 + fs1*(1.0-phase_);
        T wet = (1.0 - fade)*delay_line_.template read<I>(sd0) + fade*delay_line_.template read<I>(sd1);
        phase_ = math::wrap01(phase_ + dp1);
        p0 = math::wrap01(p0 + dp0);
        fade += fadeInc;

        if (UseInput)
        {
          binary_t is_enabled = params_.enabled.value;
          analog_t wet_level = crossfade_ = (is_enabled ? 1.0 : 0.0);
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
  T bitcrush<T, MaxBits>::process(const T& in)
  {
    rate_alpha_ += math::max<analog_t>(1.0, cast<analog_t>(rate()))*dt_;
    if (rate_alpha_ >= 1)
    {
      rate_alpha_ -= 1;
      curr_sample_ = easing::lerp(prev_input_, in, rate_alpha_);
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

  template<typename T>
  T limiter<T>::process(const T& in)
  {
    T pre = in*params_.pre_gain.value.toScale();
    T peak = math::abs(pre);
    T error = peak - params_.peak.value;
    params_.peak.value += (error > T(0)) ? T(0.05) : T(0.00002) * error;
    analog_t gain = (params_.peak.value <= 1.0 ? 1.0 : 1.0 / params_.peak.value);
    // DaisySP returns this, which sounds better for how I typically use this.
    return saturation::softlimit(pre*gain*0.7);
    // stmlib returns this, which clips more easily, but is faster.
    //return pre*gain*0.8;
  }
}
