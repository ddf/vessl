#pragma once
#include "vessl.h"

namespace vessl
{
VESSL_INLINE analog_t noise::white::operator()() const
{
  return math::random::range<analog_t>(0, 1);
}

VESSL_INLINE analog_t noise::pink::operator()()
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
      white_values_[i] = math::random::u32() % (range / count);
    }
    sum += cast<analog_t>(white_values_[i]);
  }
  max_sum_ = math::max(sum, max_sum_);
  analog_t n = sum / max_sum_;
  // vassert(!math::isNan(n) && "pink noise generated nan");
  return n;
}

VESSL_INLINE noise::red::red(analog_t sample_rate) 
: r_((sample_rate-2.0f)/sample_rate)
, rc_((1.0f - r_)*200)
, x_(0)
, y_(0)
{
  
}

VESSL_INLINE analog_t noise::red::operator()()
{
  analog_t white = math::random::range<analog_t>(-rc_, rc_);
  analog_t x = y_ + white;
  // only run the filter when we get close to going out of range
  // to compensate for wandering away from 0.
  y_ = x < -0.49 || x > 0.49f ? x - x_ + r_*y_ : x;
  x_ = x;
  return y_ + 0.5f;
}

} // namespace vessl