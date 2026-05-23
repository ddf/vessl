#pragma once

namespace vessl
{
namespace time
{
VESSL_INLINE void clockable::clock()
{
  tempo_.samples = cast<analog_t>(math::constrain(ticks_, period_min_, period_max_));
  ticks_ = 0;
  tock(0);
}

VESSL_INLINE void clockable::clock(period_t sample_delay)
{
  tempo_.samples = cast<analog_t>(math::constrain(ticks_ + sample_delay, period_min_, period_max_));
  ticks_ = 0;
  tock(sample_delay);
}
}
}