////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2019 Damien Quartz
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

#pragma once

#include <map>
#include <vector>
#include <functional>

#define _USE_MATH_DEFINES
#include <math.h>

namespace vessl
{
#pragma region interface
  template<typename T>
  class unit;

  template<typename T>
  class port
  {    
  public:
    using UT = unit<T>;

    virtual       UT& unit() = 0;
    virtual const UT& unit() const = 0;
  };

  template<typename T>
  class output : virtual public port<T>
  {
  public:
    virtual operator T() const = 0;
  };

  template<typename T>
  class input : virtual public port<T>
  {
  public:
    virtual T operator=(T v) = 0;
  };

  template<typename E>
  class array
  {
  public:
    virtual std::size_t size() const = 0;

    virtual E& operator[](std::size_t i) = 0;
    virtual const E& operator[](std::size_t i) const = 0;
  };

  template<typename T>
  class unit
  {
  public:
    using input = input<T>;
    using output = output<T>;

    virtual array<input>& inputs() = 0;
    virtual array<output>& outputs() = 0;

    // deltaTime is time in seconds since the last tick, typically 1 / sampleRate
    virtual void tick(T deltaTime) = 0;
  };

  template<typename T>
  class waveform
  {
  public:    
    /// returns the value of the waveform at phase, where phase is [0,1)
    virtual T value(T phase) const = 0;
  };

  // TODO: unpatching outputs or inputs
  template<typename T>
  class signal
  {
  public:
    using unit = unit<T>;
    using input = input<T>;
    using output = output<T>;

    signal& patch(output& src, input& dst)
    {
      //routings[&dst] = &src;

      auto currentRoute = std::find_if(routings.begin(), routings.end(), [&](auto pair) { return pair.second.dst == &dst; });
      if (currentRoute != routings.end())
      {
        routings.erase(currentRoute);
      }

      unit* u = &src.unit();
      route r = { &src, &dst };
      routings.insert(std::pair<unit*, route>(u, r));

      // add units if not already present
      if (std::find(units.begin(), units.end(), u) == units.end())
      {
        units.push_back(u);
      }

      u = &dst.unit();
      if (std::find(units.begin(), units.end(), u) == units.end())
      {
        units.push_back(u);
      }

      // sort so we get correct order for synthesis
      std::sort(units.begin(), units.end(), *this);

      return *this;
    }

    /// returns the outputs of the last unit ticked
    array<output>& tick(T deltaTime)
    {
      for (auto u : units)
      {
        u->tick(deltaTime);

        //for (auto pair : routings)
        //{
        //  if (&pair.second->unit() == u)
        //  {
        //    *pair.first = *pair.second;
        //  }
        //}

        // copy values from outputs of this unit that are patched to inputs
        auto routes = routings.equal_range(u);
        for (auto it = routes.first; it != routes.second; ++it)
        {
          *(it->second.dst) = *(it->second.src);
        }
      }

      return units.back()->outputs();
    }

    // predicate for std::sort
    bool operator()(unit* a, unit* b)
    {
      // if any output of a is patched to an input of b, a should generate before b.
      //for(auto pair : routings)
      //{
      //  if (&pair.first->unit() == b && &pair.second->unit() == a)
      //  {
      //    return true;
      //  }
      //}

      bool ab = false;
      bool ba = false;
      for (auto pair : routings)
      {
        if (pair.first == a && &pair.second.dst->unit() == b)
        {
          ab = true;
        }

        if (pair.first == b && &pair.second.dst->unit() == a)
        {
          ba = true;
        }
      }

      // only put a before b if it is patched to b and b is not patched to a in a feedback loop.
      // otherwise, our comparator will be invalid and the sort will break.
      return (ab && !ba);
    }

  private:
    std::vector<unit*> units;
    //std::map<input*, output*> routings;

    struct route
    {
      output* src;
      input* dst;
    };

    std::multimap<unit*, route> routings;
  };
#pragma endregion

#pragma region input/output  
  template<typename T>
  struct value : public input<T>, public output<T>
  {
    value() : v(0), u(0) {}
    value(UT* iu, T iv) : v(iv), u(iu) {}

    UT& unit() override { return *u; }
    const UT& unit() const override { return *u; }

    operator T() const override { return v; }
    T operator=(T i) override { v = i; return v; }

    T v;
    UT * u;
  };

  template<typename T, int N, typename E>
  struct proxyarray : public array<E>
  {
    using value = value<T>;
    proxyarray(value * const wrap) : ref(wrap) {}

    std::size_t size() const override { return N; }

    E& operator[](std::size_t i)       override { return ref[i]; }
    const E& operator[](std::size_t i) const override { return ref[i]; }

    value * const ref;
  };

  template<typename T, int I, int O>
  class io final
  {
  public:
    using input = input<T>;
    using output = output<T>;
    using value = value<T>;

    io(unit<T> * u) : in_proxy(in), out_proxy(out)
    {
      for (int i = 0; i < I; ++i)
      {
        in[i].u = u;
      }

      for (int i = 0; i < O; ++i)
      {
        out[i].u = u;
      }
    }

    value in[I];
    value out[O];

    array<input>& inputs() { return in_proxy; }
    array<output>& outputs() { return out_proxy; }

  private:
    proxyarray<T, I, input>  in_proxy;
    proxyarray<T, O, output> out_proxy;
  };
#pragma endregion

#pragma region types and utilities
  struct polarity
  {
    enum type
    {
      unipolar,
      bipolar
    };
  };

  struct interpolation
  {
    enum type
    {
      nearest,
      linear,
      cubic
    };
  };

  struct noiseTint
  {
    enum type
    {
      white,
      pink,
      red,
      brown,
    };
  };

  template<typename T>
  T clamp(T x, T a, T b)
  {
    return std::max<T>(std::min<T>(x, b), a);
  }

  template<typename T>
  T wrapPhase(T phase)
  {
    const T frac = std::modf(phase, &phase);
    return frac < 0 ? frac + 1 : frac;
  }

  template<typename T>
  T lerp(T v1, T v2, T a)
  {
    return v1 + (v2 - v1)*a;
  }

  template<typename T>
  T sample(const T* samples, T fidx, interpolation::type interp = interpolation::linear)
  {
    assert(fidx >= 0 && "fidx argument to sample must be non-negative!");

    switch (interp)
    {
      case interpolation::nearest:
        return samples[static_cast<size_t>(std::round(fidx))];

      case interpolation::linear:
      {
        T idx;
        const T frac = std::modf(fidx, &idx);
        const size_t x0 = static_cast<size_t>(idx);
        return samples[x0] + (samples[x0 + 1] - samples[x0])*frac;
      }

      case interpolation::cubic:
      {
        static const T div6 = static_cast<T>(1. / 6.);
        static const T div2 = static_cast<T>(0.5);

        T idx;
        const T f = std::modf(fidx, &idx);
        const T fm1 = f - 1.;
        const T fm2 = f - 2.;
        const T fp1 = f + 1;
        const size_t x0 = static_cast<size_t>(idx);
        return -f * fm1*fm2*div6 * samples[x0 - 1] + fp1 * fm1*fm2*div2 * samples[x0] - fp1 * f*fm2*div2 * samples[x0 + 1] + fp1 * f*fm1*div6 * samples[x0 + 2];
      }
    }

    return 0;
  }

  // a fixed-sized buffer that supports safely sampling it with a fractional index in the range [0, N-1],
  // or with a normalized "at" position in the range [0,1] where 0 will return buffer[0] and 1 will return buffer[N-1].
  template<typename T, std::size_t N, interpolation::type I = interpolation::linear>
  class samplebuffer final : public waveform<T>
  {
  public:
    samplebuffer()
    {
    }

    samplebuffer(const T(&values)[N])
    {
      for (int i = 0; i < N; ++i)
      {
        buffer[i + 1] = values[i];
      }

      // configure extra values on the ends of the buffer
      // so that sampling buffers that are periodic waveforms will work correctly
      buffer[0]   = buffer[N];
      buffer[N+1] = buffer[1];
      buffer[N+2] = buffer[2];
    }

    samplebuffer(const std::function<T(T)> generator)
    {
      T phase = 0;
      T step = (T)1 / N;
      for (int i = 0; i < N; ++i)
      {
        buffer[i + 1] = generator(phase);
        phase += step;
      }

      // configure extra values on the ends of the buffer
      // so that sampling buffers that are periodic waveforms will work correctly
      buffer[0] = buffer[N];
      buffer[N + 1] = buffer[1];
      buffer[N + 2] = buffer[2];
    }

    std::size_t size() const
    {
      return N;
    }

    T get(std::size_t i) const
    {
      return buffer[i+1];
    }

    void set(std::size_t i, T val)
    {
      buffer[i + 1] = val;
      // also update wrap-around values
      switch (i)
      {
        case 0: buffer[N+1] = buffer[1]; break;
        case 1: buffer[N+2] = buffer[2]; break;
        case N - 1: buffer[0] = buffer[N]; break;
      }
    }

    T sample(T fidx) const
    {
      return vessl::sample(buffer, fidx+1, I);
    }

    T value(T phase) const override
    {
      return sample(N*phase);
    }

  private:
    T buffer[N + 3] = {};
  };

  template<typename T, std::size_t N>
  class ringbuffer
  {
  public:
    ringbuffer()
    {
      memset(buffer, 0, sizeof(buffer));
    }

    std::size_t index(std::size_t i)
    {
      return i % N;
    }

    void push(T v)
    {
      buffer[index(write++)] = v;
    }

    T pop()
    {
      return buffer[index(read++)];
    }

    T front()
    {
      return buffer[index(read)];
    }

    void clear()
    {
      write = read = 0;
    }

    std::size_t size() const { return write - read; }

    std::size_t capacity() const { return N - size(); }

  private:
    T buffer[N];
    std::size_t write = 0;
    std::size_t read = 0;
  };  
#pragma endregion

#pragma region waveforms
  template<typename T>
  struct waves
  {
    inline static const samplebuffer<T, 1024> sine = samplebuffer<T, 1024>([](T phase) {return std::sin(phase*M_PI*2); });
  };
#pragma endregion

#pragma region unit generators
  template<typename T>
  class noise final : unit<T>
  {
  public:
    static T next(noiseTint::type tint, T dt)
    {
      switch (tint)
      {
        case noiseTint::white:
        {
          return 2 * ((T)rand() / RAND_MAX) - 1;
        }

        // #TODO: this contains clicks when run at audio frequency and I'm not sure why.
        // Need to either read up on how to do this properly and write a new implementation,
        // or find some open source code that can be used without causing license issues.
        case noiseTint::red:
        case noiseTint::brown:
        {
          static const T rc = (T)1 / (M_PI * 200);
          static const T ac = 6.2;
          static T prev = 0;
          T alpha = dt / (dt + rc);
          T white = 2 * ((T)rand() / RAND_MAX) - 1;
          prev = lerp<T>(prev, white, alpha);
          return prev * ac;
        }

        // This is the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
        // Would be good to dig into the improvements on the algorithm mentioned later in the article.
        case noiseTint::pink:
        {
          static const int range = 128;
          static const int maxKey = 0x1f;
          static int key = 0;
          static T maxSum = 90;
          static int whiteValues[6] = { rand() % (range / 6), rand() % (range / 6), rand() % (range / 6), rand() % (range / 6), rand() % (range / 6), rand() % (range / 6) };

          int lastKey = key;
          T sum = 0;
          key = key == maxKey ? 0 : ++key;

          int diff = lastKey ^ key;
          for (int i = 0; i < 6; ++i)
          {
            if ((diff & (1 << i)) != 0)
            {
              whiteValues[i] = rand() % (range / 6);
            }
            sum += whiteValues[i];
          }
          maxSum = std::max<T>(sum, maxSum);
          T n = 2 * (sum / maxSum) - 1;
          assert(!isnan(n) && "pink noise generated nan");
          return n;
        }
      }

      return 0;
    }

  public:
    noise(noiseTint::type withTint) : tint(withTint), io(this)
    {
      io.out[0] = nz[0] = nz[1] = next(withTint, 0);
      io.in[0] = 1;
    }

    noiseTint::type tint;
    input& rate = io.in[0];
    output& out = io.out[0];

    array<input>& inputs() override { return io.inputs(); }
    array<output>& outputs() override { return io.outputs(); }

    void tick(T deltaTime) override
    {
      step += deltaTime * io.in[0];
      T alpha = step / deltaTime;
      if (alpha >= 1)
      {
        nz[0] = nz[1];
        nz[1] = next(tint, deltaTime);
        alpha = wrapPhase(alpha);
        step = alpha * deltaTime;
      }
      io.out[0] = lerp(nz[0], nz[1], alpha);
    }

  private:
    T step = 0;
    T nz[2];
    io<T, 1, 1> io;
  };

  // unit that generates a linear ramp from one value to another over a duration of seconds
  template<typename T>
  class ramp final : public unit<T>
  {
  public:
    ramp(T valueBegin = 0, T valueEnd = 0, T rampDuration = 0)
      : io(this)
    {
      io.in[0] = valueBegin;
      io.in[1] = valueEnd;
      io.in[2] = rampDuration;
      io.out[0] = valueBegin;
    }

    input& begin = io.in[0];
    input& end = io.in[1];
    input& duration = io.in[2];

    output& value = io.out[0];

    array<input>& inputs() override { return io.inputs(); }
    array<output>& outputs() override { return io.outputs(); }

    bool running() const { return active; }

    void trigger()
    {
      io.out[0] = io.in[0];
      active = true;
    }

    void tick(T deltaTime)
    {
      if (active)
      {
        T dur = io.in[2];
        time += deltaTime;
        active = time < dur;
        io.out[0] = active ? lerp<T>(io.in[0], io.in[1], time / dur) : io.in[1];
      }
    }

  private:
    io<T, 3, 1> io;
    bool active;
    T time;
  };

  template<typename T, int I, int O>
  class mixer final : public unit<T>
  {
  public:
    mixer(T masterVolume)
      : audio(this)
      , volProxy(volCtrl)
      , masterCtrl(this, masterVolume)
    {
      // default all volume controls to 1
      // and the mix matrix to "unity",
      // meaning that each input is sent only to the corresponding output at full volume.
      for (int i = 0; i < I; ++i)
      {
        vol[i] = 1;
        for (int o = 0; o < O; ++o)
        {
          matrix[i][o] = T(i == o);
        }
      }
    }

    /// audio inputs
    array<input>& in = audio.inputs();
    // audio outputs
    array<output>& out = audio.outputs();

    /// mix matrix - indicates how much of each input signal to send to each output
    T matrix[I][O];

    /// volume modulation for each input
    array<input>& vol = volProxy;

    /// master volume applied to each output after accumulating inputs
    input& master = masterCtrl;

    array<input>& inputs() override { return in; }
    array<output>& outputs() override { return out; }

    void tick(T deltaTime) override
    {
      T mv = masterCtrl;
      for (int o = 0; o < O; ++o)
      {
        // initialize
        T v = 0;
        // sum all inputs to this output based on matrix and volume settings
        for (int i = 0; i < I; ++i)
        {
          v += audio.in[i] * matrix[i][o] * volCtrl[i] * mv;
        }
        audio.out[o] = v;
      }
    }

  private:
    // audio in/out
    io<T, I, O> audio;
    value<T> volCtrl[I];
    proxyarray<T, I, input> volProxy;
    value<T> masterCtrl;
  };

  // lightweight unit to transform an input value to an output via a user-provided delegate
  template<typename T, int I = 1>
  class function : public unit<T>
  {
  public:
    using delegateType = std::function<T(T, T)>;

    function(delegateType d)
      : io(this)
      , delegate(d)
    {
    }

    delegateType delegate;

    input& in = io.in[0];
    output& out = io.out[0];

    array<input>& inputs() override { return io.inputs(); }
    array<output>& outputs() override { return io.outputs(); }

    void tick(T deltaTime) override
    {
      io.out[0] = delegate(io.in[0], deltaTime);
    }

  protected:
    io<T, I, 1> io;
  };

  // multiplies 'in' by 'factor' and puts the result in 'out'
  template<typename T>
  class multiplier final : public function<T, 2>
  {
  public:
     multiplier(T amount = 0) 
      : function([this](T val, T dt){ return val * io.in[1]; })
    {
      io.in[1] = amount;
    }

    unit::input& factor = io.in[1];
  };

  // sums all inputs to a single output
  template<typename T, int I>
  class summer final : public unit<T>
  {
  public:
    summer() : io(this) { }

    array<input>& in = io.inputs();
    output& out = io.out[0];

    array<input>& inputs() override { return in; }
    array<output>& outputs() override { return io.outputs(); }

    void tick(T deltaTime) override
    {
      T sum = 0;
      for (int i = 0, sz = in.size(); i < sz; ++i)
      {
        sum += io.in[i];
      }
      io.out[0] = sum;
    }

  private:
    io<T, I, 1> io;
  };

  template<typename T>
  class oscil final : public unit<T>
  {
  public:
    using waveform = waveform<T>;

    oscil(const waveform& w, T freqInHz = 440)
      : io(this)
      , wave(w)
      , phase(0)
    {
      io.in[0] = freqInHz;
    }

    const waveform& wave;

    polarity::type polarity = polarity::bipolar;

    /// center frequency in Hz
    input& fhz = io.in[0];
    /// frequency in 1v/oct range, relative to fhz
    input& fvoct = io.in[1];
    /// current value of the oscillator
    output& out = io.out[0];

    array<input>& inputs() override { return io.inputs(); }
    array<output>& outputs() override { return io.outputs(); }

    void tick(T deltaTime) override
    {
      io.out[0] = this->polarity ? wave.value(phase) : wave.value(phase)*0.5 + 0.5;
      phase += io.in[0] * deltaTime * pow(2., io.in[1]);
      phase = wrapPhase<T>(phase);
    }

    void reset()
    {
      phase = 0;
    }

  private:
    io<T, 2, 1> io;   
    T phase;
  };

  template<typename T>
  class waveshaper final : public unit<T>
  {
  public:
    using waveform = waveform<T>;

    waveshaper(const waveform& waveShape, bool wrapInput = true)
      : shape(waveShape)
      , wrap(wrapInput)
      , io(this)
    {

    }

    // input in [-1,1] range that will be used sample the shape waveform
    input& in = io.in[0];
    // sampled shape value based on current value of in
    output& out = io.out[0];

    const waveform& shape;
    bool wrap = true;

    array<input>& inputs() override { return io.inputs(); }
    array<output>& outputs() override { return io.outputs(); }

    void tick(T deltaTime) override
    {
      assert(!isnan(io.in[0]) && "waveshaper input is nan!");
      T phase = wrap ? wrapPhase<T>(io.in[0] * 0.5 + 0.5) : clamp<T>(io.in[0] * 0.5 + 0.5, 0, 1);
      io.out[0] = shape.value(phase);
    }

  private:
    io<T, 1, 1> io;
  };

  template<typename T, std::size_t IO, std::size_t MAX_BUFFER_SIZE>
  class delay : unit<T>
  {
  public:
    delay(T timeInSeconds, T sampleRate)
      : audio(this)
      , ctrl(this)
    {
      time = timeInSeconds;
    }

    /// input signal that will be delayed
    array<input>& in = audio.inputs();
    /// output signal
    array<output>& out = audio.outputs();

    /// amount of delay time in seconds, will be clamped to MAX_BUFFER_SIZE based on current sample rate
    input& time = ctrl.in[0];
    /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
    input& feedback = ctrl.in[1];
    /// mix between dry and wet signal, clamped [0,1]
    input& mix = ctrl.in[2];

    array<input>& inputs() override { return in; }
    array<output>& outputs() override { return out; }

    void tick(T deltaTime) override
    {
      // delay time in samples, don't allow negative delay times because our buffer will fill up
      T dts = clamp<T>(ctrl.in[0] / deltaTime, 1, MAX_BUFFER_SIZE);
      T fbk = clamp<T>(ctrl.in[1], -1., 1.);
      T fade = clamp<T>(ctrl.in[2], 0., 1.);
      for (int i = 0; i < IO; ++i)
      {
        // audio input plus wet signal feedback from last tick
        T dry = audio.in[i] + ctrl.out[i] * fbk;
        if (buffer[i].capacity())
        {
          buffer[i].push(dry);
        }

        T wet = 0;

        // NOTE: this method seems to work ok for fixed delay times and slowly changing delay times,
        // but fast delay time changes generate gross artifacts,
        // which will require a resampling solution to deal with (ala VCV Rack's Delay),
        // meaning that this delay class is not entirely suitable for use in a flanger.
        //
        // how many samples do we need to pull out of the buffer to catch up with our delay time
        T read = buffer[i].size() - dts;
        if (read >= 0)
        {
          T s = buffer[i].front();
          while (read >= 1)
          {
            s = buffer[i].pop();
            read -= 1;
          }

          wet = s + (buffer[i].front() - s)*(1 - read);
        }

        ctrl.out[i] = wet;

        audio.out[i] = dry + (wet - dry)*fade;
      }
    }

  private:
    ringbuffer<T, MAX_BUFFER_SIZE> buffer[IO];
    io<T, IO, IO> audio;
    io<T, 3,  IO> ctrl; // outputs here hold the wet signal from the previous tick
  };
}

#pragma endregion

// need to break out the rotation matrix calculation into a separate class
// to allow for easier template specialization.
// but, since rotation matrices beyond 3x3 start to get real math-y,
// i'm not even sure this is something that should be included in this class.
// it might make more sense to have a separate class that is specifically about
// positioning inputs in 2d or 3d "space" parametrically,
// at let this class be a generic tool for routing inputs to outputs with scaling.
//void calcRotationMatrix()
//{
//  T rad = rotCtrl.in[0] * PI;
//  T sin = std::sin(rad);
//  T cos = std::cos(rad);
//  rotation[0][0] = cos; rotation[0][1] = -sin;
//  rotation[1][0] = sin; rotation[1][1] = cos;
//}
//
//for (int i = 0; i < I; ++i)
//{
//  T* mix = matrix[i];
//  for (int o = 0; o < O; ++o)
//  {
//    // for every in -> out mix coefficient m
//    // we need to compute m-prime,
//    // by performing a matrix multiplication across
//    // the rotation matrix's o-th row.
//    float mp = 0;
//    int ri = o;
//    for (int ro = 0; ro < O; ++ro)
//    {
//      mp += mix[ro] * rotation[ri][ro];
//    }
//
//    // we then scale the rotated mix coefficient by overall volume and apply to i-th input
//    io.out[o] = v * mp*io.in[i];
//  }
//}
