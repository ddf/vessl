#pragma once
#include <cstdint>

namespace vessl
{
  using analog_t = float;
  using digital_t = int32_t;

  // 32-bit fixed-point representation of [0,1]
  struct q32
  {
    uint32_t v_;
    constexpr q32() : v_(0u) {}
    constexpr q32(const uint32_t& v) : v_(v) {}
    constexpr q32(const q32& q) : v_(q.v_) {}
    explicit constexpr q32(const analog_t& v) : v_(0u) { *this = v; }

    constexpr q32& operator=(const q32& x) { v_ = x.v_; return *this; }
    constexpr q32& operator=(const analog_t& x) 
    { 
        analog_t f = x < -1.0f ? 1.0f : x > 1.0f ? 1.0f : x;
        if (f < 0) { f += 1.0f; }
        v_ = f * 4294967296.0;
        return *this; 
    }

    constexpr operator uint64_t() const { return v_; }

    explicit constexpr operator analog_t() const 
    {
      return v_ == UINT32_MAX ? 1.0f : v_ == 0u ? 0.f : static_cast<analog_t>(v_) / 4294967296.0f; 
    }

    // prefix increment
    q32& operator++() { v_ += 1; return *this; }
    // postfix increment
    q32 operator++(int) { q32 p = *this; operator++(); return p; } 
    // prefix decrement
    q32& operator--() { v_ -= 1; return *this; } 
    // postfix decrement
    q32 operator--(int) { q32 p = *this; operator--(); return p; }

    constexpr q32& operator+=(const q32& rhs) { v_ += rhs.v_; return *this; }
    friend constexpr q32 operator+(q32 lhs, const q32& rhs) { lhs += rhs; return lhs; }

    q32& operator-() { v_ = -v_; return *this; }
    q32& operator-=(const q32& rhs) { v_ -= rhs.v_; return *this; }
    friend q32 operator-(q32 lhs, const q32& rhs) { lhs -= rhs; return lhs; }

    q32& operator*=(const q32& rhs) 
    { 
      *this = sat(((uint64_t)v_ * (uint64_t)rhs.v_)>>31);
      return *this; 
    }
    friend q32 operator*(q32 lhs, const q32& rhs) { lhs *= rhs; return lhs; }

    friend analog_t operator*(const analog_t& lhs, const q32& rhs)
    {
      return lhs * rhs.operator vessl::analog_t();
    }

    friend analog_t operator*(const q32& lhs, const analog_t& rhs)
    {
      return lhs.operator vessl::analog_t() * rhs;
    }

    q32 scaled(int32_t factor) const
    {
      return q32(v_* factor);
    }

    q32& operator/=(const q32& rhs) 
    { 
      uint64_t n = v_; n <<= 32;
      *this = rhs.v_ == 0u ? min() : sat(n / rhs.v_); 
      return *this; 
    }
    friend q32 operator/(q32 lhs, const q32& rhs) { lhs /= rhs; return lhs; }

    friend constexpr bool operator==(const q32& lhs, const q32& rhs) { return lhs.v_ == rhs.v_; }
    friend constexpr bool operator!=(const q32& lhs, const q32& rhs) { return !(lhs == rhs); }
    friend constexpr bool operator>(const q32& lhs, const q32& rhs) { return lhs.v_ > rhs.v_; }
    friend constexpr bool operator>=(const q32& lhs, const q32& rhs) { return lhs.v_ >= rhs.v_; }
    friend constexpr bool operator<(const q32& lhs, const q32& rhs) { return lhs.v_ < rhs.v_; }
    friend constexpr bool operator<=(const q32& lhs, const q32& rhs) { return lhs.v_ <= rhs.v_; }

    static constexpr q32 max() { return UINT32_MAX; }
    static constexpr q32 mid() { return INT32_MAX; }
    static constexpr q32 min() { return 0u; }
    static constexpr q32 sat(uint64_t v) { return v > UINT32_MAX ? max() : q32((uint32_t)v); }
  };
}