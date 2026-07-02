#pragma once

namespace vessl
{
namespace sample
{

template <typename T, size_t N>
frame<T, N>::frame()
{
  as_array().fill(0);
}

template <typename T, size_t N>
frame<T, N>::frame(T m): array<T>(samples, N)
{
  as_array().fill(m);
}

template <typename T, size_t N>
frame<T, N>::frame(const frame &other): array<T>(samples, N)
{
  other.as_array().copy_to(this->as_array());
}

template <typename T, size_t N>
frame<T, N>& frame<T, N>::operator=(const frame& rhs)
{
  if (this == &rhs)
  {
    return *this;
  }
  
  rhs.as_array().copy_to(this->as_array());
  return *this;
}

template <typename T, size_t N>
frame<T, 1> frame<T, N>::to_mono() const
{
  T sum = 0;
  for (size_t c = 0; c < N; ++c)
  {
    sum += samples[c];
  }
  return frame<T, 1>(sum / N);
}

template <typename T, size_t N>
array<T> frame<T, N>::as_array() const
{
  return array<T>(samples, N);
}

template <typename T, size_t N>
matrix<T> frame<T, N>::as_matrix() const
{
  return matrix<T>(const_cast<T*>(samples), N, 1);
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator+(frame<T,N> lhs, const frame<T,N>& rhs)
{
  frame<T, N> result;
  lhs.add(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator-(frame<T,N> lhs, const frame<T,N>& rhs)
{
  frame<T, N> result;
  lhs.subtract(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(frame<T, N> lhs, const frame<T,N>& rhs)
{
  frame<T, N> result;
  lhs.multiply(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(frame<T, N> lhs, const T& rhs)
{
  frame<T, N> result;
  lhs.scale(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator*(T lhs, const frame<T, N>& rhs)
{
  frame<T, N> result;
  rhs.scale(lhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr frame<T, N> operator^(frame<T, N> lhs, const frame<T,N>& rhs)
{
  frame<T, N> result;
  for (size_t i = 0; i < N; ++i)
  {
    result[i] = cast<digital_t>(lhs[i]) ^ cast<digital_t>(rhs[i]);
  }
  return result;
}

template<typename T>
struct frame<T, 1>
{
  T samples[1];
  
  VESSL_INLINE frame() { samples[0] = 0; }
  VESSL_INLINE explicit frame(T m) { samples[0] = m; }
  VESSL_INLINE frame(const frame& other) { samples[0] = other.samples[0]; }
  VESSL_INLINE frame(frame&& other) noexcept { samples[0] = std::move(other.samples[0]); }
  VESSL_INLINE ~frame() = default;
  
  VESSL_INLINE frame& operator=(const frame& other)  // NOLINT(bugprone-unhandled-self-assignment)
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = other.samples[0];
    return *this;
  }
  
  VESSL_INLINE frame& operator=(frame&& other) noexcept
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = std::move(other.samples[0]);
    return *this;
  }
  
  VESSL_INLINE frame to_mono() const { return frame(samples[0]); }
  VESSL_INLINE array<T> as_array() const { return array<T>(samples, 1); }
  VESSL_INLINE matrix<T> as_matrix() const { return matrix<T>(samples, 1, 1); }

  VESSL_INLINE T& value() { return samples[0]; }
  VESSL_INLINE const T& value() const { return samples[0]; }
  
  VESSL_INLINE frame& operator+=(const frame& rhs)
  {
    samples[0] += rhs.samples[0];
    return *this;
  }
  
  VESSL_INLINE frame& operator-=(const frame& rhs)
  {
    samples[0] -= rhs.samples[0];
    return *this;
  }
  
  VESSL_INLINE frame& operator*=(const frame& rhs)
  {
    samples[0] *= rhs.samples[0];
    return *this;
  }
  
  VESSL_INLINE frame& operator *=(const T& rhs)
  {
    samples[0] *= rhs;
    return *this;
  }
  
  VESSL_INLINE frame& operator/=(const frame& rhs)
  {
    samples[0] /= rhs.samples[0];
    return *this;
  }
  
  VESSL_INLINE frame& operator^=(const frame& rhs)
  {
    samples[0] ^= rhs.samples[0];
    return *this;
  }
  
  VESSL_INLINE friend constexpr frame operator+(frame lhs, const frame& rhs)
  {
    return frame(lhs.value() + rhs.value());
  }
  
  VESSL_INLINE friend constexpr frame operator-(frame lhs, const frame& rhs)
  {
    return frame(lhs.value() - rhs.value());
  }
  
  VESSL_INLINE friend constexpr frame operator*(frame lhs, const frame& rhs)
  {
    return frame(lhs.value() * rhs.value());
  }
  
  VESSL_INLINE friend constexpr frame operator*(frame lhs, const T& rhs)
  {
    return frame(lhs.value() * rhs);
  }
  
  VESSL_INLINE friend constexpr frame operator*(T lhs, const frame& rhs)
  {
    return frame(lhs * rhs.value());
  }
  
  VESSL_INLINE friend constexpr frame operator^(frame lhs, const frame& rhs)
  {
    return frame(lhs.value() ^ rhs.value());
  }
};

template<typename T>
struct frame<T, 2>
{
  T samples[2];

  VESSL_INLINE frame() { samples[0] =  T(0LL); samples[1] = T(0LL); }
  VESSL_INLINE explicit frame(T m) { samples[0] = m; samples[1] = m; }
  VESSL_INLINE frame(T left, T right) { samples[0] = left, samples[1] = right; }
  VESSL_INLINE frame(const frame& other) { samples[0] = other.samples[0]; samples[1] = other.samples[1]; }
  VESSL_INLINE frame(frame&& other) noexcept { samples[0] = std::move(other.samples[0]); samples[1] = std::move(other.samples[1]); }
  VESSL_INLINE ~frame() = default;
  
  VESSL_INLINE frame& operator=(const frame& other)  // NOLINT(bugprone-unhandled-self-assignment)
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = other.samples[0];
    samples[1] = other.samples[1];
    return *this;
  }
  
  VESSL_INLINE frame& operator=(frame&& other) noexcept
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = std::move(other.samples[0]);
    samples[1] = std::move(other.samples[1]);
    return *this;
  }
  
  VESSL_INLINE frame<T, 1> to_mono() const { return frame<T, 1>((samples[0] + samples[1]) * 0.5f); }
  VESSL_INLINE array<T> as_array() const { return array<T>(samples, 2); }
  VESSL_INLINE matrix<T> as_matrix() const { return matrix<T>(samples, 2, 1); }

  VESSL_INLINE T& left() { return samples[0]; }
  VESSL_INLINE const T& left() const { return samples[0]; }
  VESSL_INLINE T& right() { return samples[1]; }
  VESSL_INLINE const T& right() const { return samples[1]; }
  
  VESSL_INLINE constexpr frame& operator+=(const frame& rhs)
  {
    *this = *this + rhs;
    return *this;
  }
  
  VESSL_INLINE constexpr frame& operator*=(const T& rhs)
  {
    samples[0] *= rhs;
    samples[1] *= rhs;
    return *this;
  }
  
  VESSL_INLINE friend constexpr frame operator+(frame lhs, const frame& rhs)
  {
    return {lhs.left() + rhs.left(), lhs.right() + rhs.right()};
  }
  
  VESSL_INLINE friend constexpr frame operator-(frame lhs, const frame& rhs)
  {
    return {lhs.left() - rhs.left(), lhs.right() - rhs.right()};
  }
  
  VESSL_INLINE friend constexpr frame operator*(frame lhs, const frame& rhs)
  {
    return {lhs.left() * rhs.left(), lhs.right() * rhs.right()};
  }
  
  VESSL_INLINE friend constexpr frame operator*(frame lhs, const T& rhs)
  {
    return {lhs.left() * rhs, lhs.right() * rhs};
  }
  
  VESSL_INLINE friend constexpr frame operator*(T lhs, const frame& rhs)
  {
    return {lhs * rhs.left(), lhs * rhs.right()};
  }
  
  VESSL_INLINE friend constexpr frame operator^(frame lhs, const frame& rhs)
  {
    return { math::xore(lhs.left(), rhs.left()), math::xore(lhs.right(), rhs.right()) };
  }
};
  
template<typename T>
struct frame<T, 3>
{
  T samples[3];

  VESSL_INLINE frame() { samples[0] =  T(0LL); samples[1] = T(0LL); samples[2] = T(0LL); }
  VESSL_INLINE explicit frame(T m) { samples[0] = m; samples[1] = m; samples[2] = m; }
  VESSL_INLINE frame(T left, T center, T right) { samples[0] = left, samples[1] = center; samples[2] = right; }
  VESSL_INLINE frame(const frame& other) { samples[0] = other.samples[0]; samples[1] = other.samples[1]; samples[2] = other.samples[2]; }
  VESSL_INLINE frame(frame&& other) noexcept { samples[0] = std::move(other.samples[0]); samples[1] = std::move(other.samples[1]); samples[2] = std::move(other.samples[2]); }
  VESSL_INLINE ~frame() = default;
    
  VESSL_INLINE frame& operator=(const frame& other)
  {
    if (this == &other)
    {
      return *this;
    }
      
    samples[0] = other.samples[0];
    samples[1] = other.samples[1];
    samples[2] = other.samples[2];
    return *this;
  }
    
  VESSL_INLINE frame& operator=(frame&& other) noexcept
  {
    if (this == &other)
    {
      return *this;
    }
      
    samples[0] = std::move(other.samples[0]);
    samples[1] = std::move(other.samples[1]);
    samples[2] = std::move(other.samples[2]);
    return *this;
  }
    
  VESSL_INLINE frame<T, 1> to_mono() const { return frame<T, 1>((samples[0] + samples[1] + samples[2]) / T(3)); }
  VESSL_INLINE array<T> as_array() const { return array<T>(const_cast<T*>(samples), 3); }
  VESSL_INLINE matrix<T> as_matrix() const { return matrix<T>(samples, 3, 1); }
  
  VESSL_INLINE T& left() { return samples[0]; }
  VESSL_INLINE const T& left() const { return samples[0]; }
  VESSL_INLINE T& center() { return samples[1]; }
  VESSL_INLINE const T& center() const { return samples[1]; }
  VESSL_INLINE T& right() { return samples[2]; }
  VESSL_INLINE const T& right() const { return samples[2]; }
  
  VESSL_INLINE T& x() { return samples[0]; }
  VESSL_INLINE const T& x() const { return samples[0]; }
  VESSL_INLINE T& y() { return samples[1]; }
  VESSL_INLINE const T& y() const { return samples[1]; }
  VESSL_INLINE T& z() { return samples[2]; }
  VESSL_INLINE const T& z() const { return samples[2]; }
  
  friend constexpr frame operator+(frame lhs, const frame& rhs)
  {
    return {lhs.left() + rhs.left(), lhs.center() + rhs.center(), lhs.right() + rhs.right()};
  }
  
  friend constexpr frame operator-(frame lhs, const frame& rhs)
  {
    return {lhs.left() - rhs.left(), lhs.center() - rhs.center(), lhs.right() - rhs.right()};
  }
  
  friend constexpr frame operator*(frame lhs, const frame& rhs)
  {
    return {lhs.left() * rhs.left(), lhs.center() * rhs.center(), lhs.right() * rhs.right()};
  }
  
  friend constexpr frame operator*(frame lhs, const T& rhs)
  {
    return {lhs.left() * rhs, lhs.center()*rhs, lhs.right() * rhs};
  }
  
  friend constexpr frame operator*(T lhs, const frame& rhs)
  {
    return {lhs * rhs.left(), lhs * rhs.center(), lhs * rhs.right()};
  }
  
  friend constexpr frame operator^(frame lhs, const frame& rhs)
  {
    return {lhs.left() ^ rhs.left(), lhs.center() ^ rhs.center(), lhs.right() ^ rhs.right()};
  }
};
} // namespace frame

VESSL_INLINE sample::gain::operator digital_t() const
{
  return cast<digital_t>(db);
}

VESSL_INLINE sample::gain::operator phase_t() const
{
  return cast<phase_t>(to_scale());
}

} // namespace vessl

template<typename T>
VESSL_INLINE T vessl::sample::interpolation::nearest::operator()(const T* buffer, analog_t frac_idx)
{
  int idx = static_cast<int>(frac_idx + 0.5f);
  return buffer[idx];
}

template<typename T>
VESSL_INLINE T vessl::sample::interpolation::linear::operator()(const T* buffer, analog_t frac_idx)
{
  int idx = static_cast<int>(frac_idx);
  T frac = cast<T>(frac_idx - idx);
  return buffer[idx] + (buffer[idx + 1] - buffer[idx]) * frac;
}

template<typename T>
VESSL_INLINE T vessl::sample::interpolation::cubic::operator()(const T* buffer, analog_t frac_idx)
{
  static constexpr analog_t div6 = (1. / 6.);
  static constexpr analog_t div2 = (0.5);

  analog_t idx;
  analog_t frc = math::mod(frac_idx, &idx);
  analog_t fm1 = frc - 1.f;
  analog_t fm2 = frc - 2.f;
  analog_t fp1 = frc + 1.f;
  size_t x0 = static_cast<size_t>(idx);
  return -frc * fm1 * fm2 * div6 * buffer[x0 - 1] 
        + fp1 * fm1 * fm2 * div2 * buffer[x0] 
        - fp1 * frc * fm2 * div2 * buffer[x0 + 1] 
        + fp1 * frc * fm1 * div6 * buffer[x0 + 2];
}

template <typename T>
VESSL_INLINE T vessl::sample::waves::bipolar::sine<T>::evaluate(phase_t phase) const
{
  return math::sin<T>(phase);
}

template <typename T>
VESSL_INLINE T vessl::sample::waves::bipolar::cosine<T>::evaluate(phase_t phase) const
{
  return math::cos<T>(phase);
}

template <typename T>
VESSL_INLINE T vessl::sample::waves::bipolar::triangle<T>::evaluate(phase_t phase) const
{
  size_t wph = static_cast<size_t>(phase) << 1;
  return wph < phase_360 ? math::lerp(T(-1), T(1), static_cast<phase_t>(wph)) 
    : math::lerp(T(1), T(-1), static_cast<phase_t>(wph - phase_360));
}

template <typename T>
T vessl::sample::waves::unipolar::sine<T>::evaluate(phase_t phase) const
{
  static constexpr T half = cast<T>(0.5f);
  return math::sin<T>(phase)*half + half;
}

template <typename T>
vessl::sample::waves::unipolar::triangle<T>::triangle()
  : attack_len_(phase_180)
  , attack_mult_(1.0f / cast<analog_t>(attack_len_))
  , decay_mult_(attack_mult_)
{
}

template <typename T>
VESSL_INLINE T vessl::sample::waves::unipolar::triangle<T>::evaluate(phase_t phase) const
{
  analog_t p = cast<analog_t>(phase);
  return phase < attack_len_ ? T(p*attack_mult_) : T((1.f - p)*decay_mult_);
}

template <typename T>
void vessl::sample::waves::unipolar::triangle<T>::set_pulse_width(phase_t pw)
{
  static constexpr phase_t pwlo = cast<phase_t>(0.01f);
  static constexpr phase_t pwhi = cast<phase_t>(0.99f);
  attack_len_ = math::constrain(pw, pwlo, pwhi);
  attack_mult_ = 1.f / cast<analog_t>(attack_len_);
  decay_mult_ = 1.f /  cast<analog_t>(phase_360 - attack_len_);
}

template <typename I, typename T>
VESSL_INLINE T vessl::sample::readf(const T *buffer, analog_t frac_idx)
{
  VASSERT(frac_idx >= 0, "fracIdx argument to sample must be non-negative");
  static I interpolator;
  return interpolator(buffer, frac_idx);
}

template <typename E, typename T>
VESSL_INLINE void vessl::sample::crossfade(const T& a, const T& b, analog_t f, T* c)
{
  static E ease;
  // copy
  *c = a;
  // in-place scale
  *c *= ease(1.0f - f);
  // constructs a T
  *c += b*ease(f);
}

template <typename T>
VESSL_INLINE void vessl::sample::crossfade(const T &a, const T &b, analog_t f, T *c)
{
  return crossfade<math::easing::linear>(a, b, f, c);
}

// tried providing this very specific specialization to see if it would improve performance 
// when using this function in Grainz on OWL2, but this doesn't appear to provide any
// substantial improvement over what gets generated by the template code.
// template <>
// VESSL_INLINE void vessl::sample::crossfade<
//   vessl::sample::frame<vessl::analog_t, 2>, 
//   vessl::math::easing::linear
// > (const frame<analog_t, 2>& a, const frame<analog_t, 2>& b, analog_t f, frame<analog_t, 2>* c)
// {
//   analog_t g = 1.f - f;
//   c->samples[0] = a.samples[0]*g + b.samples[0]*f;
//   c->samples[1] = a.samples[1]*g + b.samples[1]*f;
// }

template <typename T, vessl::size_t N, typename E>
VESSL_INLINE void vessl::sample::spatialize(const T& sample, analog_t pan, frame<T,N>* out_frame)
{
  static E ease;
  
  // see: https://www.desmos.com/calculator/atnomqw9an
  // @todo should probably just look-up standard panning/spatialization formulas.
  for (size_t i = 0; i < N; ++i)
  {
    const analog_t p = -1 + 2*static_cast<float>(i) / (N-1);
    const analog_t d = math::abs(pan - p)*0.5f;
    const analog_t a = 1.0f - d;
    out_frame->samples[i] = sample * ease(a);
  }
}

template <typename T>
VESSL_INLINE T vessl::sample::softlimit(T x)
{
  return x * (27 + x * x) / (27 + 9 * x * x);
}

template <typename T>
VESSL_INLINE T vessl::sample::softclip(T x)
{
  return x < -3 ? -1 : (x > 3 ? 1 : softlimit(x));
}

namespace vessl
{
namespace sample
{
template<typename T, size_t N, typename I>
VESSL_INLINE wavetable<T, N, I>::wavetable(source<T>& source): waveform<T>()
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
VESSL_INLINE wavetable<T, N, I>::wavetable(const waveform<T>& waveform): waveform<T>()
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
VESSL_INLINE void wavetable<T, N, I>::set(const size_t i, T val)
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
VESSL_INLINE T wavetable<T, N, I>::evaluate(phase_t phase) const
{
  analog_t idx = cast<analog_t>(phase) * N;
  return sample::readf<I>(buffer, idx);
}

template <typename T>
VESSL_INLINE ring_buffer<T>::ring_buffer(T *ring_data, size_t data_size) 
: array<T>(ring_data, data_size)
, write_index_(0)
{
  
}

template<typename T>
VESSL_INLINE T ring_buffer<T>::write(const T& v)
{
  T o = data_[write_index_];
  data_[write_index_++] = v;
  if (write_index_ == size_)
  {
    write_index_ = 0;
  }
  return o;
}

template <typename T>
VESSL_INLINE size_t ring_buffer<T>::get_write_index() const
{
  return write_index_;
}

template <typename T>
VESSL_INLINE void ring_buffer<T>::set_write_index(size_t index)
{
  write_index_ = index%size_;
}

template<typename T>
VESSL_INLINE ring_buffer<T> ring_buffer<T>::operator<<(typename array<T>::reader r)
{
  VASSERT(r.available() < size(), "reader size is larger than ring size");
  while (r)
  {
    write(r.read());
  }
  return *this;
}

template <typename T>
VESSL_INLINE T ring_buffer<T>::overdub(const T &v, size_t write_offset)
{
  size_t idx = get_write_index() + write_offset;
  if (idx >= size_) idx %= size_;
  T& s = data_[idx];
  s += v;
  return s;
}

template <typename T>
void ring_buffer<T>::overwrite(const T &v, size_t write_offset)
{
  size_t idx = get_write_index() + write_offset;
  if (idx >= size_) idx %= size_;
  data_[idx] = v;
}

template<typename T>
VESSL_INLINE T delay_line<T>::read(size_t sample_delay) const
{
  //VASSERT(sample_delay < size());
  size_t sz = size();
  size_t idx = get_write_index() - sample_delay;
  return idx >= sz ? data()[idx%sz] : data()[idx];
}

template<typename T>
VESSL_INLINE T delay_line<T>::readf(analog_t sample_delay) const
{
  const size_t sz = size();
  const analog_t widx = get_write_index();
  const analog_t idx = sample_delay <= widx ? widx - sample_delay : widx + sz - sample_delay;
  const size_t lidx = static_cast<size_t>(idx);
  const size_t hidx = (lidx+1) % sz;
  const analog_t t = idx - lidx;
  return math::lerp(data()[lidx], data()[hidx], t);
}

template<typename T>
VESSL_INLINE T delay_line<T>::evaluate(phase_t phase) const
{
  analog_t size_f = cast<analog_t>(size());
  analog_t sample_delay = cast<analog_t>(phase_360 - phase) * size_f;
  return readf(sample_delay);
}

} // namesapce sample
} // namespace vessl