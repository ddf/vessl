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

#ifdef round
#undef round
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

#include "vessl_qmath.h"

// Note: In all classes using typename T, it is assumed to be POD and to have support for all arithmetic operators
namespace vessl
{
using char_t    = char;
using size_t    = uint64_t;
using binary_t  = bool;
using digital_t = int64_t;
using analog_t  = float;
using phase_t   = uint32_t;
using q31_t     = q31;

static constexpr phase_t phase_360  = UINT32_MAX;
static constexpr phase_t phase_180  = UINT32_MAX >> 1;
static constexpr phase_t phase_90   = phase_180 >> 1;
static constexpr phase_t phase_270  = phase_180 + phase_90;
static constexpr phase_t phase_zero = 0UL;

// we use this in place of static_cast throughout the library for non-pointer types
// so that we can specialize conversions between some of our value types (e.g. phase_t <--> analog_t)
template<typename T, typename F>
VESSL_INLINE constexpr T cast(const F& from) { return static_cast<T>(from); }

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
VESSL_INLINE T mod(T v, T* i) { return ::modf(v, i); }

template<typename T>
VESSL_INLINE T pow(T x, T y) { return ::pow(x, y); }

template<typename T>
VESSL_INLINE T round(T x) { return ::round(x); }
      
template<typename T>
VESSL_INLINE T floor(T x) { return ::floor(x); }

template<typename T, typename R = T>
VESSL_INLINE T sin(R r) { return ::sin(r); }

template<typename T, typename R = T>
VESSL_INLINE T cos(R r) { return ::cos(r); }

template<typename T>
VESSL_INLINE T sqrt(T x) { return ::sqrt(x); }
      
template<typename T>
VESSL_INLINE T sqrt2() { static T v = sqrt(2); return v; }

template<typename T>
VESSL_INLINE T tan(T x) { return ::tan(x); }
    
analog_t decibels_to_scale(analog_t db);

analog_t scale_to_decibels(analog_t scale);

template<typename T, typename D>
T lerp(T begin, T end, D t);

template<typename T>
T wrap(T val, T low, T high);

template<typename T>
VESSL_INLINE T wrap01(T val) { return wrap(val, T(0), T(1)); }

template<typename T>
VESSL_INLINE binary_t is_nan(T n) { return isnan(n); }
      
// xore because xor is a keyword
template<typename T>
VESSL_INLINE T xore(const T& a, const T& b) { return a ^ b; }
  
namespace random
{
// implements the xorshifter algorithm
static constexpr uint32_t u32_max = UINT32_MAX;
static uint32_t ru32_seed = 33641;
VESSL_INLINE void su32(uint32_t seed) {ru32_seed = seed; }
uint32_t u32();

template<typename T>
T range(T low, T high);
} // namespace random
  
namespace easing
{
struct linear
{
  template<typename D>
  D operator()(D t) const;
};

struct smoothstep
{
  analog_t operator()(analog_t t) const;
};

namespace quad
{
struct in { analog_t operator()(analog_t t) const; };
struct out { analog_t operator()(analog_t t) const; };
struct in_out { analog_t operator()(analog_t t) const; };
struct out_in { analog_t operator()(analog_t t) const; };
}

namespace expo
{
struct in { analog_t operator()(analog_t t) const; };
struct out { analog_t operator()(analog_t t) const; }; 
struct in_out { analog_t operator()(analog_t t) const; };
}

template<typename T>
T smooth(T value, T target, analog_t degree = 0.9f);
    
template<typename T>
struct smoother
{
  T        value;
  analog_t degree;
      
  explicit smoother(analog_t smoothing_degree = 0.9f, T initial_value = T(0));

  VESSL_INLINE explicit operator T() const { return value; }
      
  // so we can use this like OWL's SmoothValue
  smoother& operator=(const T& v);
};
} // namespace easing

// easing is first parameter so T can be deduced
template<typename E, typename T, typename D>
T interp(T begin, T end, D t);

} // namespace math

// analog unipolar noise generators that generate values in the range [0,1]
namespace noise
{
struct white
{
  explicit white(analog_t sample_rate) { (void)sample_rate; }
  analog_t operator()() const;
};

// Implements the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
// Would be good to dig into the improvements on the algorithm mentioned later in the article.
struct pink
{
  explicit pink(analog_t sample_rate) { (void)sample_rate; }
      
  analog_t operator()();

private:
  static constexpr int range = 128;
  static constexpr int count = 6;
  static constexpr int max_key = 0x1f;

  int key_ = 0;
  analog_t max_sum_ = 90;
  uint32_t white_values_[count] = {
    math::random::u32() % (range / count), math::random::u32() % (range / count), math::random::u32() % (range / count),
    math::random::u32() % (range / count), math::random::u32() % (range / count), math::random::u32() % (range / count)
  };
};

// Brownian noise (i.e. random wander) run thru a DC blocking filter.
// See: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
// @todo still a bit crunchy
struct red
{
  explicit red(analog_t sample_rate);

  analog_t operator()();

private:
  analog_t r_;
  analog_t rc_;
  analog_t x_, y_;
};
} // namespace noise

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

  [[nodiscard]] virtual binary_t is_full() const = 0;
  virtual void write(const T& value) = 0;

  VESSL_INLINE sink& operator<<(const T& value) { write(value); return *this; }
  VESSL_INLINE explicit operator binary_t() const { return !is_full(); }
};

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
    VESSL_INLINE explicit reader() 
    : source<T>()
    , begin_(nullptr), head_(nullptr), end_(nullptr) {}
      
    VESSL_INLINE reader(const T* data, size_t size) 
    : source<T>()
    , begin_(data), head_(data), end_(data + size) {}
      
    VESSL_INLINE explicit reader(array source) : reader(source.data_, source.size_) {}

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
  void copy_to(array dest) const;

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

template<typename T>
class matrix;

namespace sample
{
// a struct to hold one sample frame
template<typename T, size_t N>
struct frame : array<T>
{
  T samples[N];

  frame();
  explicit frame(T m);
  frame(const frame& other);
  frame(frame&&) = default;
  frame& operator=(const frame&) = default;
  frame& operator=(frame&&) = default;
  ~frame() = default;

  // mixdown to a mono frame
  frame<T, 1> to_mono() const;
          
  // matrix view of this frame
  matrix<T> as_matrix() const;
};

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator+(frame<T,N> lhs, const frame<T,N>& rhs);

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator-(frame<T,N> lhs, const frame<T,N>& rhs);

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(frame<T, N> lhs, const frame<T,N>& rhs);

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(frame<T, N> lhs, const T& rhs);

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(T lhs, const frame<T, N>& rhs);

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator^(frame<T, N> lhs, const frame<T,N>& rhs);
    
template<typename T>
struct type
{
  using mono = frame<T,1>;
  using stereo = frame<T,2>;
  using lcr = frame<T,3>;
  using vector3 = frame<T,3>;
};
  
// stored as decibels (0 = unity gain).
struct gain
{
  analog_t db;
          
  VESSL_INLINE gain() : db(0) {}
  VESSL_INLINE explicit gain(analog_t dbv) : db(dbv) {}
  VESSL_INLINE static gain from_scale(analog_t scale) { return gain(math::scale_to_decibels(scale)); }
  VESSL_INLINE static gain from_decibels(analog_t db) { return gain(db); }

  // implement casting operators so that gain can be used as a parameter type.
  VESSL_INLINE explicit operator binary_t() const  { return db >= 0; }
  VESSL_INLINE explicit operator digital_t() const;
  VESSL_INLINE explicit operator analog_t() const  { return db; }
  VESSL_INLINE explicit operator phase_t() const;

  VESSL_INLINE analog_t to_scale() const {  return math::decibels_to_scale(db); }
  VESSL_INLINE analog_t to_decibels() const { return db; }
};
    
// interface for a waveform that can be evaluated using a normalized phase value
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
  VESSL_INLINE T evaluate(phase_t phase) const override;
};

template<typename T>
struct cosine final : waveform<T>
{
  VESSL_INLINE T evaluate(phase_t phase) const override;
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
}
  
template<typename T, typename I = interpolation::linear<T>>
T read_interpolated(const T* buffer, analog_t frac_idx);
  
template<typename T, typename E = math::easing::linear>
T crossfade(T a, T b, analog_t f);
  
// lovingly borrowed from pichenettes/stmlib
template<typename T>
T softlimit(T x);

// lovingly borrowed from pichenettes/stmlib
template<typename T>
T softclip(T x);
    
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

// aka circular array
template<typename T>
class ring_buffer : array<T>
{
public:
  ring_buffer(T* ring_data, size_t data_size) : array<T>(ring_data, data_size), head_(ring_data + data_size - 1) {}

  // expose direct access to underlying array data
  using array<T>::data;
  using array<T>::size;

  void write(const T& v);
  size_t get_write_index() const { return head_ - array<T>::data_; }
  void set_write_index(size_t index) { head_ = array<T>::data_ + index%array<T>::size_; }

  ring_buffer operator<<(typename array<T>::reader r);
      
private:
  T* head_;
};
  
template<typename T>
class delay_line final : public ring_buffer<T>, public waveform<T>
{
public:
  delay_line(T* delay_line_data, size_t data_size) : ring_buffer<T>(delay_line_data, data_size) {}

  using ring_buffer<T>::data;
  using ring_buffer<T>::size;
  using ring_buffer<T>::get_write_index;
  using ring_buffer<T>::set_write_index;

  // reads behind the write head with sampleDelay (i.e. the ith sample previously written)
  // where a delay of 0 samples will give the most recently written value.
  T read(size_t sample_delay) const;

  // reads behind the write head with a fractional sampleDelay and given interpolation
  template<typename I>
  T read(analog_t sample_delay) const;

  // phase will be wrapped to [-1,1] where 0 is the oldest sample recorded
  T evaluate(phase_t phase) const override;
};
} // namespace sample

namespace filtering
{
using gain_t = sample::gain;

// common filter Q values for biquad filters
namespace q
{
template<typename T>
VESSL_INLINE constexpr T butterworth() { return cast<T>(0.70710678118); } // 1/sqrt(2)

template<typename T>
VESSL_INLINE constexpr T sallen_key() { return cast<T>(0.5); } 

template<typename T>
VESSL_INLINE constexpr T bessel() { return cast<T>(0.57735026919); } // 1/sqrt(3)
}
  
struct args
{
  analog_t sr;
  analog_t hz;
  analog_t q;
  gain_t   g;
      
  VESSL_INLINE args(analog_t sample_rate, analog_t hertz, analog_t kyu, gain_t gain) 
  : sr(sample_rate), hz(hertz), q(kyu), g(gain)
  {
  }
      
  // helper for biquad
  VESSL_INLINE analog_t omega() const { return hz * math::pi<analog_t>() / sr; }
};
    
// DC blocking filter, see: https://www.dsprelated.com/freebooks/filters/DC_Blocker.html
template<typename T>
struct dc_block
{
  T x1 = T(0), y1 = T(0);
  void process(const T* source, T* dest, size_t block_size, const args& args);
};
    
template<typename T>
struct data
{
  using coeff_t = analog_t;
  using state_t = T;
      
  coeff_t* coeff;
  state_t* state;

  data(coeff_t* coeff_data, size_t coeff_size, T* state_data, size_t state_size);
};

// based on https://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
template<size_t Stages>
struct biquad
{
  static constexpr size_t coeff_num = 5;
      
  template<typename T, size_t States>
  struct cascade : data<T>
  {
    static constexpr size_t coeff_count = coeff_num*Stages;
    static constexpr size_t state_count = States*Stages;
        
    cascade() : data<T>(co, coeff_count, st, state_count), co{} {}
        
    analog_t co[coeff_count];
    T st[state_count];
  };
      
  template<typename T, class CoGen>
  // ReSharper disable once CppInconsistentNaming
  struct df2t final : cascade<T, 2>
  {
    static CoGen cg;
        
    using cascade<T,2>::coeff_count;
    using cascade<T,2>::state_count;
         
    void process(const T* source, T* dest, size_t block_size, const args& args);
    // ReSharper disable once CppMemberFunctionMayBeStatic
    [[nodiscard]] size_t stage_count() const { return Stages; }
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
} // namespace filtering
  

// matrix_data is a separate struct so that we can specialize for ARM easily
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
    
  VESSL_INLINE void clear();
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

// not even sure this belongs in here
template<typename T>
class transform33 : public matrix<T>
{
public:
  transform33() : matrix<T>(data_, 3, 3)
  {
    set_identity();
  }

  transform33(const transform33& other) : matrix<T>(data_, 3, 3)
  {
    memcpy(data_, other.data_, sizeof(T)*3*3);
  }
    
  VESSL_INLINE void set_identity() 
  {
    matrix<T>::clear();
    for (size_t i = 0; i < 3; i++) {
      matrix<T>::set(i,i, T(1LL));
    }
  }
      
  void set_euler(phase_t pitch, phase_t yaw, phase_t roll);

  VESSL_INLINE void set_euler_radians(analog_t pitch_radians, analog_t yaw_radians, analog_t roll_radians)
  {
    return set_euler(cast<phase_t>(pitch_radians / math::two_pi<analog_t>()), 
      cast<phase_t>(yaw_radians / math::two_pi<analog_t>()), 
      cast<phase_t>(roll_radians / math::two_pi<analog_t>()));
  }
      
  [[nodiscard]] VESSL_INLINE sample::frame<T,3> multiply(const sample::frame<T,3>& input)
  {
    using m = matrix<T>;
    sample::frame<T,3> output;

    output[0] = m::get(0,0) * input[0] + m::get(0,1) * input[1] + m::get(0,2) * input[2];
    output[1] = m::get(1,0) * input[0] + m::get(1,1) * input[1] + m::get(1,2) * input[2];
    output[2] = m::get(2,0) * input[0] + m::get(2,1) * input[1] + m::get(2,2) * input[2];

    // this might be faster?
    //mtrx.multiply(input.toMatrix(), output.toMatrix());

    return output;
  }
      
private:
  T data_[3 * 3];
};
  
template<typename T>
class list
{
public:
  list() = default;
  virtual ~list() = default;
  list(const list&) = default;
  list(list&&) = default;
  list& operator=(const list&) = default;
  list& operator=(list&&) = default;
    
  [[nodiscard]] virtual size_t size() const = 0;
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
VESSL_INLINE typename list<T>::iterator begin(const list<T>& lst) { return list<T>::iterator::begin(lst); }
  
template<typename T>
VESSL_INLINE typename list<T>::iterator end(const list<T>& lst) { return list<T>::iterator::end(lst); }

namespace time
{
// can be used by units as an indication for how to treat changes in time values (see: delay & freeze)
enum class mode : uint8_t
{
  snap, // use duration value directly
  slew, // smooth duration when it changes
  fade, // crossfade to new duration across a block of samples
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

  static duration from_bpm(analog_t bpm, analog_t sample_rate) { return duration(sample_rate/(bpm*b_to_f)); }
  static duration from_seconds(analog_t seconds, analog_t sample_rate) { return duration(sample_rate*seconds); }
  [[nodiscard]] analog_t to_bpm(analog_t sample_rate) const { return f_to_b*(sample_rate/samples); }
  [[nodiscard]] analog_t to_seconds(analog_t sample_rate) const { return samples/sample_rate; }
};
  
  
// classes can subclass this to add support for tempo detection of a clock signal (i.e. pulse train)
class clockable
{
public:
  using period_t = uint32_t;

  clockable(analog_t sample_rate, period_t sample_period_min, period_t sample_period_max, analog_t bpm = 60)
  : tempo_(duration::from_bpm(bpm, sample_rate)), period_min_(sample_period_min), period_max_(sample_period_max)
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
      
  duration tempo_;
  period_t period_min_;
  period_t period_max_;
  period_t ticks_;
  analog_t sample_rate_;
};
} // namespace time

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
      
    // note to self: additional built-ins should not be structs or classes
    // doing so means that if a user of the library wants to use an enum as a parameter,
    // they will need to specialize the cast function to provide a conversion from all
    // struct/class built-ins to their enum type, otherwise they will get a confusing compilation error.

    // user provided type, stored as a void*, must be convertible to all other parameter types.
    // even if the conversion is meaningless.
    user = UINT8_MAX
  };
    
  template<typename T>
  static constexpr value_type type_of() { return value_type::user; }
    
  typedef uint32_t id_t;
    
  struct desc
  {
    const char_t* name;
    id_t          id;
    value_type    type;
      
    VESSL_INLINE constexpr desc(const char_t* n, id_t i, value_type t) : name(n), id(i), type(t) {}
    VESSL_INLINE static constexpr desc empty() { return {"", 0, value_type::none}; }
  };
    
  template<size_t N>
  struct desc_list
  {
    static constexpr size_t size = N;
    desc list[N];
    constexpr desc operator[](id_t id) const;
  };
    
  template<typename T>
  struct data
  {
    T value = T(0);
  };
    
  template<typename T>
  constexpr parameter(const desc& param_desc, const data<T>& param_data) : desc_(param_desc), data_(&const_cast<data<T>&>(param_data).value) {}
      
  // explicitly declared copy-constructors so that they will be used instead of the copy-assignment override
  // when returning parameter objects by value (as in plist::element_at implementations)
  constexpr parameter(parameter& param) : desc_(param.desc_), data_(param.data_) {};
  constexpr parameter(const parameter& param) : desc_(param.desc_), data_(const_cast<void*>(param.data_)) {};
    
  [[nodiscard]] VESSL_INLINE constexpr const desc& description() const { return desc_; }
    
  template<typename T>
  T read() const;

  VESSL_INLINE binary_t   read_binary()  const { return read<binary_t>(); }
  VESSL_INLINE digital_t  read_digital() const { return read<digital_t>(); }
  VESSL_INLINE analog_t   read_analog()  const { return read<analog_t>(); }
  VESSL_INLINE phase_t    read_phase()   const { return read<phase_t>(); }
    
  VESSL_INLINE explicit operator binary_t()  const { return read<binary_t>(); }
  VESSL_INLINE explicit operator digital_t() const { return read<digital_t>(); }
  VESSL_INLINE explicit operator analog_t()  const { return read<analog_t>(); }
  VESSL_INLINE explicit operator phase_t()   const { return read<phase_t>(); }

  // static-cast T to parameter type before assign
  template<typename T>
  parameter& write(const T& value);

  template<typename T>
  VESSL_INLINE parameter& operator=(const T& value) { write(value); return *this; }
    
  parameter& operator=(const parameter& rhs);

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
  
template<typename T>
struct param : parameter::data<T>
{
  constexpr parameter operator()(const char_t* name, parameter::id_t id) const;
};

using gain_t = sample::gain;
using duration_t = time::duration;
  
typedef param<analog_t>    analog_p;
typedef param<digital_t>   digital_p;
typedef param<binary_t>    binary_p;
typedef param<phase_t>     phase_p;
typedef param<gain_t>      gain_p;
typedef param<duration_t>  duration_p;

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
  [[nodiscard]] virtual desc description() const { return { "", nullptr, 0 }; }
  [[nodiscard]] virtual const parameter_list& parameters() const = 0;
    
  // @todo get parameter by name / id
};
  
VESSL_INLINE const parameter::desc* begin(const unit::desc& desc) { return desc.params; }
VESSL_INLINE const parameter::desc* end(const unit::desc& desc) { return desc.params + desc.param_count; }

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
    
  struct sample
  {
    processor* proc;
    I input;
      
    sample(processor& processor, I sample) : proc(&processor), input(sample) {}
    I& operator>>(O& out) { out = proc->process(input); return out; }
  };
    
  struct block
  {
    processor* proc;
    array<I> input;
      
    block(processor& processor, array<I> buffer) : proc(&processor), input(buffer) {}
    // implies in-place processing of block (but what if it didn't? how do?)
    block operator>>(processor& next) { proc->process(input, input); return block(next, input); }
    array<I> operator>>(array<O> out) { proc->process(input, out); return out; }
  };
    
  struct stream : source<O>
  {
    processor* proc;
    source<I>* src;
      
    stream(processor& processor, source<I>& input) : proc(&processor), src(&input) {}
    [[nodiscard]] binary_t is_empty() const override { return src->is_empty(); }
    O read() override { return proc->process(src->read()); }
    O& operator>>(O& rhs) { rhs = read(); return rhs; }
    stream operator>>(processor& processor) { return stream(processor, *this); }
    sink<I>& operator>>(sink<O>& out) { proc->process(*src, out); return out; }
    array<I> operator>>(array<O> out)
    {
      auto w = out.make_writer(); proc->process(*src, w); return out;
    }
  };
};

template<typename I, typename O>
typename processor<I,O>::sample operator>>(const parameter& p, processor<I,O>& proc)
{
  return typename processor<I,O>::sample(proc, p.read<I>());
}

template<typename I, typename O>
typename processor<I,O>::block operator>>(array<I> in, processor<I,O>& proc)
{
  return typename processor<I,O>::block(proc, in);
}

template<typename I, typename O>
typename processor<I,O>::block operator>>(const I& in, processor<I,O>& proc)
{
  return typename processor<I,O>::block(proc, sample::frame<I,1>(in));
}

template<typename I, typename O>
typename processor<I,O>::stream operator>>(source<I>& in, processor<I,O>& proc)
{
  return typename processor<I,O>::stream(proc, in);
}

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
} // namespace vessl

#include "vessl_core.inl"
#include "vessl_math.inl"
#include "vessl_qmath.inl"
#include "vessl_noise.inl"
#include "vessl_sample.inl"
#include "vessl_filtering.inl"
#include "vessl_time.inl"
#include "vessl_parameter.inl"

#include "vessl_units.h"