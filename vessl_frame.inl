#pragma once

namespace vessl
{
namespace frame
{

template <typename T, size_t N>
channels<T, 1> channels<T, N>::toMono() const
{
  T sum = 0;
  for (size_t c = 0; c < N; ++c)
  {
    sum += samples[c];
  }
  return channels<T, 1>(sum / N);
}

template <typename T, size_t N>
matrix<T> channels<T, N>::toMatrix() const
{
  return matrix<T>(const_cast<T*>(samples), N, 1);
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator+(channels<T,N> lhs, const channels<T,N>& rhs)
{
  channels<T, N> result;
  lhs.add(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator-(channels<T,N> lhs, const channels<T,N>& rhs)
{
  channels<T, N> result;
  lhs.subtract(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator*(channels<T, N> lhs, const channels<T,N>& rhs)
{
  channels<T, N> result;
  lhs.multiply(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator*(channels<T, N> lhs, const T& rhs)
{
  channels<T, N> result;
  lhs.scale(rhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator*(T lhs, const channels<T, N>& rhs)
{
  channels<T, N> result;
  rhs.scale(lhs, result);
  return result;
}

template<typename T, size_t N>
VESSL_INLINE constexpr channels<T, N> operator^(channels<T, N> lhs, const channels<T,N>& rhs)
{
  channels<T, N> result;
  for (size_t i = 0; i < N; ++i)
  {
    result[i] = cast<digital_t>(lhs[i]) ^ cast<digital_t>(rhs[i]);
  }
  return result;
}

template<typename T>
struct channels<T, 1> : array<T>
{
  T samples[1];
  
  VESSL_INLINE channels() : array<T>(samples, 1) { samples[0] = 0; }
  VESSL_INLINE explicit channels(T m) : array<T>(samples, 1) { samples[0] = m; }
  VESSL_INLINE channels(const channels& other) : array<T>(samples, 1) { samples[0] = other.samples[0]; }
  VESSL_INLINE channels(channels&& other) noexcept : array<T>(samples, 1) { samples[0] = std::move(other.samples[0]); }
  VESSL_INLINE ~channels() = default;
  
  VESSL_INLINE channels& operator=(const channels& other)  // NOLINT(bugprone-unhandled-self-assignment)
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = other.samples[0];
    return *this;
  }
  
  VESSL_INLINE channels& operator=(channels&& other) noexcept
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = std::move(other.samples[0]);
    return *this;
  }
  
  VESSL_INLINE channels toMono() const { return channels(samples[0]); }
  VESSL_INLINE matrix<T> toMatrix() const { return matrix<T>(samples, 1, 1); }

  VESSL_INLINE T& value() { return samples[0]; }
  VESSL_INLINE const T& value() const { return samples[0]; }
};

template<typename T>
struct channels<T, 2> : array<T>
{
  T samples[2];

  VESSL_INLINE channels() : array<T>(samples, 2) { samples[0] =  T(0LL); samples[1] = T(0LL); }
  VESSL_INLINE explicit channels(T m) : array<T>(samples, 2) { samples[0] = m; samples[1] = m; }
  VESSL_INLINE channels(T left, T right) : array<T>(samples, 2) { samples[0] = left, samples[1] = right; }
  VESSL_INLINE channels(const channels& other) : array<T>(samples, 2) { samples[0] = other.samples[0]; samples[1] = other.samples[1]; }
  VESSL_INLINE channels(channels&& other) noexcept : array<T>(samples, 2) { samples[0] = std::move(other.samples[0]); samples[1] = std::move(other.samples[1]); }
  VESSL_INLINE ~channels() = default;
  
  VESSL_INLINE channels& operator=(const channels& other)  // NOLINT(bugprone-unhandled-self-assignment)
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = other.samples[0];
    samples[1] = other.samples[1];
    return *this;
  }
  
  VESSL_INLINE channels& operator=(channels&& other) noexcept
  {
    if (this == &other)
    {
      return *this;
    }
    
    samples[0] = std::move(other.samples[0]);
    samples[1] = std::move(other.samples[1]);
    return *this;
  }
  
  VESSL_INLINE channels<T, 1> toMono() const { return channels<T, 1>((samples[0] + samples[1]) * 0.5f); }
  VESSL_INLINE matrix<T> toMatrix() const { return matrix<T>(samples, 2, 1); }

  VESSL_INLINE T& left() { return samples[0]; }
  VESSL_INLINE const T& left() const { return samples[0]; }
  VESSL_INLINE T& right() { return samples[1]; }
  VESSL_INLINE const T& right() const { return samples[1]; }
};
  
template<typename T>
struct channels<T, 3> : array<T>
{
  T samples[3];

  VESSL_INLINE channels() : array<T>(samples, 3) { samples[0] =  T(0LL); samples[1] = T(0LL); samples[2] = T(0LL); }
  VESSL_INLINE explicit channels(T m) : array<T>(samples, 3) { samples[0] = m; samples[1] = m; samples[2] = m; }
  VESSL_INLINE channels(T left, T center, T right) : array<T>(samples, 2) { samples[0] = left, samples[1] = center; samples[2] = right; }
  VESSL_INLINE channels(const channels& other) : array<T>(samples, 2) { samples[0] = other.samples[0]; samples[1] = other.samples[1]; samples[2] = other.samples[2]; }
  VESSL_INLINE channels(channels&& other) noexcept : array<T>(samples, 2) { samples[0] = std::move(other.samples[0]); samples[1] = std::move(other.samples[1]); samples[2] = std::move(other.samples[2]); }
  VESSL_INLINE ~channels() = default;
    
  VESSL_INLINE channels& operator=(const channels& other)  // NOLINT(bugprone-unhandled-self-assignment)
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
    
  VESSL_INLINE channels& operator=(channels&& other) noexcept
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
    
  VESSL_INLINE channels<T, 1> toMono() const { return channels<T, 1>((samples[0] + samples[1] + samples[2]) / T(3)); }
  VESSL_INLINE matrix<T> toMatrix() const { return matrix<T>(samples, 3, 1); }
  
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
};

} // namespace frame
} // namespace vessl
