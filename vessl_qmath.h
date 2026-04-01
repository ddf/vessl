/*
MIT License

Copyright (c) 2024 Vaclav Mach (Bastl Instruments)
Copyright (c) 2024 Marek Mach (Bastl Instruments)
Copyright (c) 2026 Damien Quartz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <cstdint>

// Fixed-point types in struct form based on the qmath types in Bastl Instrument's Kastle2 framework.
namespace vessl
{
  using analog_t = float;
  
  // 32-bit fixed-point type
  struct q31
  {
    int32_t v_;
    VESSL_INLINE constexpr q31() : v_(0) {}
    VESSL_INLINE constexpr q31(const q31& q) = default;
    VESSL_INLINE explicit constexpr q31(const int32_t& v) : v_(v) {}

    VESSL_INLINE constexpr q31& operator=(const q31& x) = default;

    // prefix increment
    VESSL_INLINE q31& operator++() { v_ += 1; return *this; }
    // postfix increment
    VESSL_INLINE q31 operator++(int) { q31 p = *this; operator++(); return p; } 
    // prefix decrement
    VESSL_INLINE q31& operator--() { v_ -= 1; return *this; } 
    // postfix decrement
    VESSL_INLINE q31 operator--(int) { q31 p = *this; operator--(); return p; }

    VESSL_INLINE constexpr q31& operator+=(const q31& rhs) 
    { 
      v_ = sat(static_cast<int64_t>(v_) + static_cast<int64_t>(rhs.v_)).v_;
      return *this; 
    }
    VESSL_INLINE friend constexpr q31 operator+(q31 lhs, const q31& rhs) { lhs += rhs; return lhs; }

    VESSL_INLINE q31& operator-() { v_ = (v_ == INT32_MIN ? INT32_MAX : -v_); return *this; }

    VESSL_INLINE q31& operator-=(const q31& rhs) 
    {
      v_ = sat(static_cast<int64_t>(v_) - static_cast<int64_t>(rhs.v_)).v_; 
      return *this; 
    }
    VESSL_INLINE friend q31 operator-(q31 lhs, const q31& rhs) { lhs -= rhs; return lhs; }

    VESSL_INLINE q31& operator*=(const q31& rhs) 
    { 
      v_ = sat((static_cast<int64_t>(v_) * static_cast<int64_t>(rhs.v_))>>31).v_;
      return *this; 
    }
    VESSL_INLINE friend q31 operator*(q31 lhs, const q31& rhs) { lhs *= rhs; return lhs; }

    [[nodiscard]] VESSL_INLINE q31 scaled(int32_t factor) const
    {
      return sat(static_cast<int64_t>(v_) * factor);
    }

    [[nodiscard]] VESSL_INLINE q31 scaled(analog_t factor) const
    {
      return sat(static_cast<int64_t>(static_cast<analog_t>(v_) * factor));
    }

    VESSL_INLINE q31& operator/=(const q31& rhs) 
    { 
      int64_t n = static_cast<int64_t>(v_) << 31;
      v_ = (rhs.v_ == 0L ? mid() : sat(n / rhs.v_)).v_; 
      return *this; 
    }
    VESSL_INLINE friend q31 operator/(q31 lhs, const q31& rhs) { lhs /= rhs; return lhs; }

    VESSL_INLINE friend constexpr bool operator==(const q31& lhs, const q31& rhs) { return lhs.v_ == rhs.v_; }
    VESSL_INLINE friend constexpr bool operator!=(const q31& lhs, const q31& rhs) { return !(lhs == rhs); }
    VESSL_INLINE friend constexpr bool operator>(const q31& lhs, const q31& rhs) { return lhs.v_ > rhs.v_; }
    VESSL_INLINE friend constexpr bool operator>=(const q31& lhs, const q31& rhs) { return lhs.v_ >= rhs.v_; }
    VESSL_INLINE friend constexpr bool operator<(const q31& lhs, const q31& rhs) { return lhs.v_ < rhs.v_; }
    VESSL_INLINE friend constexpr bool operator<=(const q31& lhs, const q31& rhs) { return lhs.v_ <= rhs.v_; }

    VESSL_INLINE static constexpr q31 max() { return q31(INT32_MAX); }
    VESSL_INLINE static constexpr q31 half() { return q31(0x3fffffffL); }
    VESSL_INLINE static constexpr q31 mid() { return {}; }
    VESSL_INLINE static constexpr q31 min() { return q31(INT32_MIN); }
    VESSL_INLINE static constexpr q31 sat(int64_t v) { return v > INT32_MAX ? max() : v < INT32_MIN ? min() : q31(static_cast<int32_t>(v)); }
    VESSL_INLINE static constexpr q31 recip(const analog_t a) { return a == 0 ? mid() : sat(INT32_MAX / a); }
  };
}