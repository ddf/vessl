#pragma once
#include <cstdint>

namespace vessl
{
  using analog_t = float;

  // 32-bit fixed-point representation of [-1,1]
  struct q31
  {
    int32_t v_;
    constexpr q31() : v_(0) {}
    constexpr q31(const q31& q) : v_(q.v_) {}
    explicit constexpr q31(const int32_t& v) : v_(v) {}
    explicit constexpr q31(const int64_t& v) : v_(0) { *this = sat(v); }
    explicit constexpr q31(const analog_t& v) : v_(0) { *this = v; }

    constexpr q31& operator=(const q31& x) { v_ = x.v_; return *this; }
    constexpr q31& operator=(const analog_t& x) 
    { 
        analog_t f = x < -1.0f ? -1.0f : x > 1.0f ? 1.0f : x;
        v_ = static_cast<int32_t>(f * 2147483648.0f);
        return *this;
    }

    explicit constexpr operator analog_t() const 
    {
      return static_cast<analog_t>(v_) / 2147483648.0f; 
    }

    // prefix increment
    q31& operator++() { v_ += 1; return *this; }
    // postfix increment
    q31 operator++(int) { q31 p = *this; operator++(); return p; } 
    // prefix decrement
    q31& operator--() { v_ -= 1; return *this; } 
    // postfix decrement
    q31 operator--(int) { q31 p = *this; operator--(); return p; }

    constexpr q31& operator+=(const q31& rhs) 
    { 
      v_ = sat((int64_t)v_ + (int64_t)rhs.v_).v_;
      return *this; 
    }
    friend constexpr q31 operator+(q31 lhs, const q31& rhs) { lhs += rhs; return lhs; }

    // adds this and rhs, allowing the value to wrap within our min/max range.
    // useful when using q31 as phase, as we do throughout the library
    q31& spill(const q31& rhs)
    {
      int64_t s = (int64_t)v_ + (int64_t)rhs.v_;
      v_ = s > INT32_MAX ? (int32_t)(s - INT32_MAX) : s < INT32_MIN ? (int32_t)(s + INT32_MAX) : s;
      return *this;
    }

    static constexpr q31 spill(q31 lhs, const q31& rhs)
    {
      return lhs.spill(rhs);
    }

    q31& operator-() { v_ = (v_ == INT32_MIN ? INT32_MAX : INT32_MAX - v_); return *this; }
    q31& operator-=(const q31& rhs) 
    {
      v_ = sat((int64_t)v_ - (int64_t)rhs.v_).v_; 
      return *this; 
    }
    friend q31 operator-(q31 lhs, const q31& rhs) { lhs -= rhs; return lhs; }

    q31& operator*=(const q31& rhs) 
    { 
      v_ = sat(((int64_t)v_ * (int64_t)rhs.v_)>>31).v_;
      return *this; 
    }
    friend q31 operator*(q31 lhs, const q31& rhs) { lhs *= rhs; return lhs; }

    friend analog_t operator*(const analog_t& lhs, const q31& rhs)
    {
      return lhs * rhs.operator vessl::analog_t();
    }

    friend analog_t operator*(const q31& lhs, const analog_t& rhs)
    {
      return lhs.operator vessl::analog_t() * rhs;
    }

    q31 scaled(int32_t factor) const
    {
      return q31(v_* factor);
    }

    q31& operator/=(const q31& rhs) 
    { 
      int64_t n = (int64_t)v_ << 31;
      v_ = (rhs.v_ == (int32_t)0u ? mid() : sat(n / rhs.v_)).v_; 
      return *this; 
    }
    friend q31 operator/(q31 lhs, const q31& rhs) { lhs /= rhs; return lhs; }

    friend constexpr bool operator==(const q31& lhs, const q31& rhs) { return lhs.v_ == rhs.v_; }
    friend constexpr bool operator!=(const q31& lhs, const q31& rhs) { return !(lhs == rhs); }
    friend constexpr bool operator>(const q31& lhs, const q31& rhs) { return lhs.v_ > rhs.v_; }
    friend constexpr bool operator>=(const q31& lhs, const q31& rhs) { return lhs.v_ >= rhs.v_; }
    friend constexpr bool operator<(const q31& lhs, const q31& rhs) { return lhs.v_ < rhs.v_; }
    friend constexpr bool operator<=(const q31& lhs, const q31& rhs) { return lhs.v_ <= rhs.v_; }

    static constexpr q31 max() { return q31(INT32_MAX); }
    static constexpr q31 mid() { return q31(); }
    static constexpr q31 min() { return q31(INT32_MIN); }
    static constexpr q31 sat(int64_t v) { return v > INT32_MAX ? max() : v < INT32_MIN ? min() : q31((int32_t)v); }
    static constexpr q31 recip(const analog_t a) { return a == 0 ? mid() : sat((int64_t)INT32_MAX / a); }
  };
}