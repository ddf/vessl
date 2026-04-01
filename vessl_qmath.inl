/*
MIT License

Copyright (c) 2024 Vaclav Mach (Bastl Instruments)
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

#include "vessl.h"
#include "vessl_qmath.h"
#include "vessl_qmath_lut.h"

// followed by template specializations for the vessl:q31 type.
namespace vessl
{
  template<>
  [[nodiscard]] VESSL_INLINE constexpr q31 cast<q31, analog_t>(const analog_t& from)
  {
    analog_t f = from < -1.0f ? -1.0f : from > 1.0f ? 1.0f : from;
    int32_t v = static_cast<int32_t>(f * 2147483648.0f);
    return q31(v);
  }
  
  template<>
  [[nodiscard]] VESSL_INLINE constexpr analog_t cast<analog_t, q31>(const q31& from)
  {
    return static_cast<analog_t>(from.v_) / 2147483648.0f;
  }

  template<>
  [[nodiscard]] VESSL_INLINE constexpr q31 cast<q31, phase_t>(const phase_t& from)
  {
    return q31(static_cast<int32_t>(from / 2));
  }

  template<>
  [[nodiscard]] VESSL_INLINE constexpr phase_t cast<phase_t, q31>(const q31& from)
  {
    int32_t v = from.v_ < 0 ? from.v_ + INT32_MAX + 1 : from.v_;
    // off by one when v == INT32_MAX, but probably OK?
    return static_cast<uint32_t>(v) << 1;
  }

  [[nodiscard]] VESSL_INLINE constexpr analog_t operator*(const analog_t& lhs, const q31& rhs)
  {
    return lhs * cast<analog_t>(rhs);
  }

  [[nodiscard]] VESSL_INLINE constexpr analog_t operator*(const q31& lhs, const analog_t& rhs)
  {
    return cast<analog_t>(lhs) * rhs;
  }
    
  namespace math
  {
    template<>
    [[nodiscard]] VESSL_INLINE q31 abs<q31>(const q31& x)
    {
      if (x == q31::min())
      {
          return q31::max();
      }
      return x.v_ > 0 ? x : q31(-x.v_);
    }

    template<>
    [[nodiscard]] VESSL_INLINE q31 sin<q31, phase_t>(phase_t phase) 
    { 
      int32_t v = static_cast<int32_t>(phase >> 1);
      return q31(qmath_sine_table[v >> QMATH_SINE_TABLE_SHIFT_Q31]);
    }

    template<>
    [[nodiscard]] VESSL_INLINE q31 cos<q31, phase_t>(phase_t phase) 
    { 
      return sin<q31>(phase + PHASE_90);
    }

    template<>
    [[nodiscard]] VESSL_INLINE q31 sin<q31, analog_t>(analog_t radians) 
    { 
      return sin<q31>(cast<phase_t>(radians / twoPi<analog_t>()));
    }

    template<>
    [[nodiscard]] VESSL_INLINE q31 cos<q31, analog_t>(analog_t radians) 
    { 
      return cos<q31>(cast<phase_t>(radians / twoPi<analog_t>()));
    }
  }

  namespace easing
  {
    template<>
    [[nodiscard]] VESSL_INLINE q31 lerpp<q31>(q31 begin, q31 end, phase_t t)
    {
      return t == PHASE_ZERO ? begin 
          : t == PHASE_360 ? end 
          : begin + (end-begin)*cast<q31>(t);
    }
  }
}
