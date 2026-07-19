#pragma once

namespace vessl
{
namespace time
{
VESSL_INLINE duration::operator bool() const
{
  return math::abs(samples) >= math::epsilon<analog_t>();
}

VESSL_INLINE clockable::clockable(analog_t sample_rate, period_t sample_period_min, period_t sample_period_max, analog_t bpm)
  : tempo_(duration::from_bpm(bpm, sample_rate))
  , period_min_(sample_period_min)
  , period_max_(sample_period_max)
  , ticks_(sample_period_max) // start unclocked
  , sample_rate_(sample_rate)
{}

VESSL_INLINE void clockable::clock()
{
  tempo_.samples = static_cast<analog_t>(math::constrain(ticks_, period_min_, period_max_));
  ticks_ = 0;
  tock(0);
}

VESSL_INLINE void clockable::clock(period_t sample_delay)
{
  tempo_.samples = static_cast<analog_t>(math::constrain(ticks_ + sample_delay, period_min_, period_max_));
  ticks_ = 0;
  tock(sample_delay);
}

VESSL_INLINE bool clockable::is_clocked() const
{
  return ticks_ < period_max_;
}

VESSL_INLINE void clockable::tick()
{
  ticks_ = math::min(ticks_ + 1,  period_max);
}

VESSL_INLINE void clockable::tick(period_t t)
{
  ticks_ = ticks_ < period_max - t ? ticks_ + t : period_max;
}

VESSL_INLINE void clockable::tock(period_t sample_delay)
{
  (void)sample_delay;
}
}
}