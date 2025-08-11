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
#include <cassert>

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace vessl
{
  template<typename T>
  class array
  {
  protected:
    T* data;
    size_t size;
  public:
    array() : data(nullptr), size(0) {}
    array(T* data, size_t size) : data(data), size(size) {}

    T* getData() { return data; }
    const T* getData() const { return data; }
    size_t getSize() const { return size; }
    bool isEmpty() const { return size == 0; }
    T& operator[](size_t index) { return data[index]; }
    const T& operator[](size_t index) const { return data[index]; }

    // ranged-based for support
    T* begin() { return data; }
    T* end() { return data + size; }

    class reader
    {
    protected:
      const T* begin;
      const T* head;
      const T* end;
    public:
      reader(T* data, size_t size) : begin(data), head(data), end(data + size) {}

      size_t available() const { return end - head; }
      explicit operator bool() const { return head != end; }

      const T& peek() const { return *head; }
      const T& operator*() const { return *head; }
      
      const T& read() { return *head++; }

      reader& reset() { head = begin; return *this; }
    };

    class writer
    {
    protected:
      T* head;
      const T* end;
    public:
      writer(T* data, size_t size) : head(data), end(data + size) {}

      size_t available() const { return end - head; }
      void write(const T& v) { *head++ = v; }

      template<typename R>
      writer operator<<(R& r)
      {
        size_t sz = r.available();
        assert(available() >= sz);
        size_t blocks = sz >> 2u;
        sz -= blocks*4;
        T a, b, c, d;
        while (blocks--)
        {
          a = r.read(); b = r.read(); c = r.read(); d = r.read();
          write(a); write(b); write(c); write(d);
        }
        while (sz--)
        {
          write(r.read());
        }
        return writer(head, available());
      }
    };
  };

  template<typename T>
  typename array<T>::writer& operator<<(typename array<T>::writer& w, const T& v)
  {
    w.write(v);
    return w;
  }

  template<typename T>
  class ring : array<T>, array<T>::writer, array<T>::reader
  {
  public:
    ring(T* inData, size_t inSize) : array<T>(inData, inSize), array<T>::writer(inData, inSize), array<T>::reader(inData, inSize)
    {
    }

    size_t getWriteCapacity() const
    {
      return array<T>::writer::head > array<T>::reader::head
      ? array<T>::writer::available() + (array<T>::reader::head - array<T>::data)
      : array<T>::reader::head - array<T>::writer::head;
    }
    
    void write(const T& v)
    {
      assert(getWriteCapacity() > 0);
      array<T>::writer::write(v);
      // reset the writer when we get to the end of our array
      if (array<T>::writer::available() == 0)
      {
        array<T>::writer::head = array<T>::data;
      }
    }

    ring& operator<<(typename array<T>::reader& r)
    {
      assert(r.available() < getWriteCapacity());
      while (r)
      {
        write(r.read());
      }
      return *this;
    }

    size_t getReadCapacity() const
    {
      return array<T>::writer::head >= array<T>::reader::head
      ? array<T>::writer::head - array<T>::reader::head
      : array<T>::reader::available() + (array<T>::writer::head - array<T>::data);
    }

    // fulfill the reader contract for block copy
    size_t available() const { return getReadCapacity(); }

    const T& read()
    {
      assert(getReadCapacity() > 0);
      const T& result = array<T>::reader::read();
      if (array<T>::reader::available() == 0)
      {
        array<T>::reader::reset();
      }
      return result;
    }
  };
  
  template<typename T>
  class parameter
  {
  public:
    parameter() : n(""), v(0) {}
    parameter(const char * name, T value) : n(name), v(value) {}
    
    const char* name() const { return n; }

    // for the sane
    T read() const { return v; }
    parameter& write(T value) { v = value; return *this; }

    // just how much silly syntactic sugar can we add here, lol
    const T& operator*() const { return v; }
    operator T*() { return &v; }
    
    T operator!() const { return v == 0 ? T(1) : T(1) / v; }
    T operator-() const { return -v; }

  private:
    const char* n;
    T v;
  };

  template<>
  inline float parameter<float>::operator!() const { return v > 0 && v < FLT_EPSILON ? 1.0f : 1.0f / v; }

  template<>
  inline double parameter<double>::operator!() const { return v > 0 && v < DBL_EPSILON ? 1.0 : 1.0 / v; }

  template<>
  inline long double parameter<long double>::operator!() const { return v > 0 && v < LDBL_EPSILON ? 1.0 : 1.0 / v; }

  template<typename T>
  parameter<T>& operator<<(parameter<T>& p, const T& v)
  {
    p.write(v);
    return p;
  }

  template<typename T>
  parameter<T>& operator>>(parameter<T>& p, T& v)
  {
    v = p.read();
    return p;
  }

  template<typename T>
  class generator
  {
  public:
    generator() = default;
    virtual ~generator() = default;
    generator(const generator&) = delete;
    generator(const generator&&) = delete;
    generator& operator=(const generator&) = delete;
    generator& operator=(generator&&) = delete;
    
    virtual T generate() = 0;  // NOLINT(portability-template-virtual-member-function)
  };

  template<typename T>
  class processor
  {
  public:
    processor() = default;
    virtual ~processor() = default;
    processor(const processor&) = delete;
    processor(const processor&&) = delete;
    processor& operator=(const processor&) = delete;
    processor& operator=(processor&&) = delete;
    
    virtual void process(const T& in, T* out) = 0;
  };

  template<typename T, size_t N>
  struct uinit
  {
    const char* name;
    parameter<T> params[N];
  };
  
  template<typename T>
  class unit : array<parameter<T>>
  {
  public:
    template<size_t N>
    explicit unit(uinit<T, N>& init, float sampleRate = 1) : array<parameter<T>>(init.params, N), n(init.name) { setSampleRate(sampleRate); }
    unit(const char* name, parameter<T>* params, int paramsCount, float sampleRate = 1) : array<parameter<T>>(params, paramsCount), n(name) { setSampleRate(sampleRate); }
    virtual ~unit() = default;
    unit(const unit&) = default;
    unit(unit&&) = default;
    unit& operator=(const unit&) = default;
    unit& operator=(unit&&) = default;
    
    const char* name() const { return n; }
    
    float getSampleRate() const { return sampleRate; }
    
    void setSampleRate(float sr)
    {
      sampleRate = sr;
      deltaTime = 1.0f/sr;
      onSampleRateChanged();
    }

    parameter<T>& getParameter(size_t index) { return this->operator[](index); }
    const parameter<T>& getParameter(size_t index) const { return this->operator[](index); }
    size_t getParameterCount() const { return array<parameter<T>>::getSize(); }

    using array<parameter<T>>::begin;
    using array<parameter<T>>::end;

  protected:
    virtual void onSampleRateChanged() {}
    float dt() const { return deltaTime; }
    
  private:
    const char* n;
    float sampleRate;
    float deltaTime;
  };

  template<typename T>
  class unitGenerator : public unit<T>, public generator<T>
  {
  public:
    template<size_t N>
    explicit unitGenerator(uinit<T, N>& init, float sampleRate = 1) : unit<T>(init, sampleRate), generator<T>() {}
    unitGenerator(const char* name, parameter<T>* params, int paramsCount, float sampleRate = 1) : unit<T>(name, params, paramsCount, sampleRate), generator<T>() {}
  };

  template<typename T>
  class unitProcessor : public unit<T>, public processor<T>
  {
  public:
    template<size_t N>
    explicit unitProcessor(uinit<T, N>& init, float sampleRate = 1) : unit<T>(init, sampleRate), processor<T>() {}
    unitProcessor(const char* name, parameter<T>* params, int paramsCount, float sampleRate = 1) : unit<T>(name, params, paramsCount, sampleRate), processor<T>() {}
  };
  
  namespace interpolation
  {
    template<typename T>
    struct nearest
    {
      T operator()(const T* buffer, double fracIdx)
      {
        return buffer[static_cast<size_t>(roundf(fracIdx))];
      }
    };

    template<typename T>
    struct linear
    {
      T operator()(const T* buffer, double fracIdx)
      {
        T idx;
        T frac = modf(fracIdx, &idx);
        size_t x0 = static_cast<size_t>(idx);
        return buffer[x0] + (buffer[x0 + 1] - buffer[x0])*frac;
      }
    };

    template<typename T>
    struct cubic
    {
      T operator()(const T* buffer, double fracIdx)
      {
        static const T DIV6 = static_cast<T>(1. / 6.);
        static const T DIV2 = static_cast<T>(0.5);

        double idx;
        double f = modf(fracIdx, &idx);
        double fm1 = f - 1.;
        double fm2 = f - 2.;
        double fp1 = f + 1;
        size_t x0 = static_cast<size_t>(idx);
        return -f * fm1*fm2*DIV6 * buffer[x0 - 1] + fp1 * fm1*fm2*DIV2 * buffer[x0] - fp1 * f*fm2*DIV2 * buffer[x0 + 1] + fp1 * f*fm1*DIV6 * buffer[x0 + 2];
      }
    };
  };
  
  template<typename T, typename I = interpolation::linear<T>>
  T sample(const T* buffer, double fracIdx)
  {
    assert(fracIdx >= 0 && "fracIdx argument to sample must be non-negative!");
    static I interpolator;
    return interpolator(buffer, fracIdx);
  }

  template<typename T>
  T wrap(T val, T low, T high)
  {
    T diff = high - low;
    while (val < low)
    {
      val += diff;
    }
    while (val > high)
    {
      val -= diff;
    }
    return val;
  }

  template<typename T>
  T wrap01(T val) { return wrap(val, T(0), T(1)); }

  template<>
  inline float wrap01(float v) { float i; return modf(v, &i); }

  // a waveform that can be evaluated using a phase value, which will be wrapped to the range [0,1)
  template<typename T>
  class waveform
  {
  public:
    waveform() = default;
    virtual ~waveform() = default;
    waveform(const waveform&) = delete;
    waveform(const waveform&&) = delete;
    waveform& operator=(const waveform&) = delete;
    waveform& operator=(waveform&&) = delete;
    
    T operator()(double phase) const { return eval(wrap01(phase)); }
    
  protected:
    virtual T eval(double phase) const = 0;  // NOLINT(portability-template-virtual-member-function)
  };
  
  namespace waves
  {
    template<typename T = double>
    class sine final : public waveform<T>
    {
    public:
      sine() = default;
    protected:
      T eval(double phase) const override { return sin(2*PI*phase); }
    };
    
    template<>
    inline float sine<float>::eval(double phase) const { return sinf(static_cast<float>(2 * PI * phase)); }
  }
  
  // a fixed-sized buffer that supports safely sampling it with a fractional index in the range [0, N-1],
  // or with a normalized position (i.e. phase) in the range [0,1] where 0 will return buffer[0] and 1 will return buffer[N-1].
  template<typename T, size_t N, typename I = interpolation::linear<T>>
  class wavetable final : public waveform<T>
  {
  public:
    wavetable(): waveform<T>() {}

    wavetable(const T(&values)[N]): waveform<T>()
    {
      for (int i = 0; i < N; ++i)
      {
        buffer[i + 1] = values[i];
      }

      // configure extra values on the ends of the buffer
      // so that sampling buffers that are periodic waveforms will work correctly
      buffer[0] = buffer[N];
      buffer[N + 1] = buffer[1];
      buffer[N + 2] = buffer[2];
    }

    wavetable(const waveform<T>& generator): waveform<T>()
    {
      double phase = 0;
      double step = T(1) / N;
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

    // ReSharper disable once CppMemberFunctionMayBeStatic
    size_t size() const
    {
      return N;
    }

    T get(const size_t i) const
    {
      return buffer[i+1];
    }

    void set(const size_t i, T val)
    {
      buffer[i + 1] = val;
      // also update wrap-around values
      switch (i)
      {
        case 0: buffer[N+1] = buffer[1]; break;
        case 1: buffer[N+2] = buffer[2]; break;
        case N - 1: buffer[0] = buffer[N]; break;
        default: break;
      }
    }

    T sample(double fracIdx) const
    {
      return vessl::sample(buffer, fracIdx+1);
    }

  protected:
    // implement waveform
    T eval(double phase) const override { return sample(phase*N); }

  private:
    T buffer[N + 3] = {};
  };

  // @todo easings namespace similar to interpolation namespace
  template<typename T>
  T lerp(T v1, T v2, T a)
  {
    return v1 + (v2 - v1)*a;
  }
  
  enum class noiseTint : uint8_t
  {
    white,
    pink,
    red,
    brown,
  };
  
  template<typename T>
  class noise final : public unitGenerator<T>
  {
    uinit<T, 2> init = {
      "noise", {{ "tint", 0 }, { "rate", 1 }}
    };

    using unit<T>::dt;
    
  public:
    explicit noise(float sampleRate, noiseTint withTint = noiseTint::white)
    : unitGenerator<T>(init, sampleRate), step(0)
    {
      nz[0] = nz[1] = next(withTint, 0);
      tint().write(static_cast<float>(withTint));
    }
    noise(const noise&) = default;
    noise(noise&&) = default;
    noise& operator=(const noise&) = default;
    noise& operator=(noise&&) = default;
    ~noise() override = default;

    parameter<T>& tint() { return init.params[0]; }
    parameter<T>& rate() { return init.params[1]; }

    static T next(noiseTint tint, T dt);

    T generate() override
    {
      step += dt() * rate().read();
      float alpha = step / dt();
      if (alpha >= 1)
      {
        noiseTint tnt = static_cast<noiseTint>(tint().read());
        nz[0] = nz[1];
        nz[1] = next(tnt, dt());
        alpha = wrap01(alpha);
        step = alpha * dt();
      }
      return lerp(nz[0], nz[1], alpha);
    }

  private:
    float step;
    T nz[2];
  };

  template<typename T>
  T noise<T>::next(noiseTint tint, T dt) {
    switch (tint)
    {
      case noiseTint::white:
      {
        return 2 * (static_cast<T>(rand()) / RAND_MAX) - 1;
      }

      // #TODO: this contains clicks when run at audio frequency and I'm not sure why.
      // Need to either read up on how to do this properly and write a new implementation,
      // or find some open source code that can be used without causing license issues.
      case noiseTint::red:
      case noiseTint::brown:
      {
        static const T RC = static_cast<T>(1) / (PI * 200);
        static const T AC = 6.2;
        static T prev = 0;
        T alpha = dt / (dt + RC);
        T white = 2 * (static_cast<T>(rand()) / RAND_MAX) - 1;
        prev = lerp<T>(prev, white, alpha);
        return prev * AC;
      }

      // This is the Voss algorithm (see: http://www.firstpr.com.au/dsp/pink-noise/)
      // Would be good to dig into the improvements on the algorithm mentioned later in the article.
      case noiseTint::pink:
      {
        static constexpr int RANGE = 128;
        static constexpr int MAX_KEY = 0x1f;
        static int key = 0;
        static T maxSum = 90;
        static int whiteValues[6] = { rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6), rand() % (RANGE / 6) };

        int lastKey = key;
        T sum = 0;
        key = key == MAX_KEY ? 0 : ++key;

        int diff = lastKey ^ key;
        for (int i = 0; i < 6; ++i)
        {
          if ((diff & (1 << i)) != 0)
          {
            whiteValues[i] = rand() % (RANGE / 6);
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

  
  // unit that generates a linear ramp from one value to another over a duration of seconds
  // @todo implement easings above and add that as a template parameter
  template<typename T>
  class ramp final : public unitGenerator<T>
  {
    uinit<T, 4> init {
      "ramp", {{"begin", 0}, {"end", 0}, {"duration", 0}, {"eor", 0}}
    };

    using unit<T>::dt;
    
  public:
    explicit ramp(float sampleRate, float durationInSeconds = 0, float fromValue = 0, float toValue = 0) : unitGenerator<T>(init, sampleRate)
    , t(0)
    {
      // this is fun
      from() << fromValue;
      to() << toValue;
      duration() << durationInSeconds;
    }
    ramp(const ramp&) = default;
    ramp(ramp&&) = default;
    ramp& operator=(const ramp&) = default;
    ramp& operator=(ramp&&) = default;
    ~ramp() override = default;

    // ins
    parameter<T>& from() { return init.params[0]; }
    parameter<T>& to() { return init.params[1]; }
    parameter<T>& duration() { return init.params[2]; }

    // outs
    const parameter<T>& eor() const { return init.params[3]; }
    // could also add t as an out.
      
    bool isActive() const { return *eor() < 1; }
  
    void trigger()
    {
      if (*duration() > 0)
      {
        t = 0;
        eor_() << T(0);
      }
      else
      {
        t = 1;
        eor_() << T(1);
      }
    }
  
    T generate() override
    {
      float lt = t;
      if (isActive())
      {
        t += dt() * !duration();
        if (t >= 1)
        {
          eor_() << T(1);
        }
      }
      return lerp<T>(*from(), *to(), lt);
    }
  
  private:
    parameter<T>& eor_() { return init.params[3]; }
    float t;
  };
  //
  // template<typename T, int I, int O>
  // class mixer final : public unit<T>
  // {
  // public:
  //   mixer(T masterVolume)
  //     : audio(this)
  //     , volProxy(volCtrl)
  //     , masterCtrl(this, masterVolume)
  //   {
  //     // default all volume controls to 1
  //     // and the mix matrix to "unity",
  //     // meaning that each input is sent only to the corresponding output at full volume.
  //     for (int i = 0; i < I; ++i)
  //     {
  //       vol[i] = 1;
  //       for (int o = 0; o < O; ++o)
  //       {
  //         matrix[i][o] = T(i == o);
  //       }
  //     }
  //   }
  //
  //   /// audio inputs
  //   array<input>& in = audio.inputs();
  //   // audio outputs
  //   array<output>& out = audio.outputs();
  //
  //   /// mix matrix - indicates how much of each input signal to send to each output
  //   T matrix[I][O];
  //
  //   /// volume modulation for each input
  //   array<input>& vol = volProxy;
  //
  //   /// master volume applied to each output after accumulating inputs
  //   input& master = masterCtrl;
  //
  //   array<input>& inputs() override { return in; }
  //   array<output>& outputs() override { return out; }
  //
  //   void tick(T deltaTime) override
  //   {
  //     T mv = masterCtrl;
  //     for (int o = 0; o < O; ++o)
  //     {
  //       // initialize
  //       T v = 0;
  //       // sum all inputs to this output based on matrix and volume settings
  //       for (int i = 0; i < I; ++i)
  //       {
  //         v += audio.in[i] * matrix[i][o] * volCtrl[i] * mv;
  //       }
  //       audio.out[o] = v;
  //     }
  //   }
  //
  // private:
  //   // audio in/out
  //   io<T, I, O> audio;
  //   value<T> volCtrl[I];
  //   proxyarray<T, I, input> volProxy;
  //   value<T> masterCtrl;
  // };
  //
  // // lightweight unit to transform an input value to an output via a user-provided delegate
  // template<typename T, int I = 1>
  // class function : public unit<T>
  // {
  // public:
  //   using delegateType = std::function<T(T, T)>;
  //
  //   function(delegateType d)
  //     : io(this)
  //     , delegate(d)
  //   {
  //   }
  //
  //   delegateType delegate;
  //
  //   input& in = io.in[0];
  //   output& out = io.out[0];
  //
  //   array<input>& inputs() override { return io.inputs(); }
  //   array<output>& outputs() override { return io.outputs(); }
  //
  //   void tick(T deltaTime) override
  //   {
  //     io.out[0] = delegate(io.in[0], deltaTime);
  //   }
  //
  // protected:
  //   io<T, I, 1> io;
  // };
  //
  // // multiplies 'in' by 'factor' and puts the result in 'out'
  // template<typename T>
  // class multiplier final : public function<T, 2>
  // {
  // public:
  //    multiplier(T amount = 0) 
  //     : function([this](T val, T dt){ return val * io.in[1]; })
  //   {
  //     io.in[1] = amount;
  //   }
  //
  //   unit::input& factor = io.in[1];
  // };
  //
  // // sums all inputs to a single output
  // template<typename T, int I>
  // class summer final : public unit<T>
  // {
  // public:
  //   summer() : io(this) { }
  //
  //   array<input>& in = io.inputs();
  //   output& out = io.out[0];
  //
  //   array<input>& inputs() override { return in; }
  //   array<output>& outputs() override { return io.outputs(); }
  //
  //   void tick(T deltaTime) override
  //   {
  //     T sum = 0;
  //     for (int i = 0, sz = in.size(); i < sz; ++i)
  //     {
  //       sum += io.in[i];
  //     }
  //     io.out[0] = sum;
  //   }
  //
  // private:
  //   io<T, I, 1> io;
  // };

  template<typename T, typename W>
  class oscil final : public unitGenerator<T>
  {
    W wave;
    T phase;
    T dt;

    uinit<T, 4> init = {
      "oscil", {{ "frequency", 0 }, { "fm (lin)", 0 }, { "fm (v/oct)", 0 }, { "phase mod", 0 }}
    };
    
  public:
    explicit oscil(T sampleRate, T freqInHz = 440) : unitGenerator<T>(init)
    , phase(0) , dt(1.0f/sampleRate)
    {
      fHz().write(freqInHz);
    }
    oscil(const oscil&) = default;
    oscil(oscil&&) = default;
    oscil& operator=(const oscil&) = default;
    oscil& operator=(oscil&&) = default;
    ~oscil() override = default;

    W& waveform() { return wave; }
    const W& waveform() const { return wave; }

    // frequency in Hz without FM applied
    parameter<T>& fHz() { return init.params[0]; }
    // linear frequency modulation
    parameter<T>& fmLin() { return init.params[1]; }
    // v/oct (exponential) frequency modulation
    parameter<T>& fmExp() { return init.params[2]; }
    // phase modulation
    parameter<T>& pm() { return init.params[3]; }
    
    T generate() override
    {
      T val = wave(phase + pm().read()); 
      phase += (fHz().read()*exp2(fmExp().read()) + fmLin().read())*dt;
      phase = wrap01(phase);
      return val;
    }

    void reset()
    {
      phase = 0;
    }
  };
  
//
//   template<typename T>
//   class waveshaper final : public unit<T>
//   {
//   public:
//     using waveform = waveform<T>;
//
//     waveshaper(const waveform& waveShape, bool wrapInput = true)
//       : shape(waveShape)
//       , wrap(wrapInput)
//       , io(this)
//     {
//
//     }
//
//     // input in [-1,1] range that will be used sample the shape waveform
//     input& in = io.in[0];
//     // sampled shape value based on current value of in
//     output& out = io.out[0];
//
//     const waveform& shape;
//     bool wrap = true;
//
//     array<input>& inputs() override { return io.inputs(); }
//     array<output>& outputs() override { return io.outputs(); }
//
//     void tick(T deltaTime) override
//     {
//       assert(!isnan(io.in[0]) && "waveshaper input is nan!");
//       T phase = wrap ? wrapPhase<T>(io.in[0] * 0.5 + 0.5) : clamp<T>(io.in[0] * 0.5 + 0.5, 0, 1);
//       io.out[0] = shape.value(phase);
//     }
//
//   private:
//     io<T, 1, 1> io;
//   };
//
//   template<typename T, std::size_t IO, std::size_t MAX_BUFFER_SIZE>
//   class delay : unit<T>
//   {
//   public:
//     delay(T timeInSeconds, T sampleRate)
//       : audio(this)
//       , ctrl(this)
//     {
//       time = timeInSeconds;
//     }
//
//     /// input signal that will be delayed
//     array<input>& in = audio.inputs();
//     /// output signal
//     array<output>& out = audio.outputs();
//
//     /// amount of delay time in seconds, will be clamped to MAX_BUFFER_SIZE based on current sample rate
//     input& time = ctrl.in[0];
//     /// amount of signal to feedback, can be negative to invert feedback signal, clamped [-1,1]
//     input& feedback = ctrl.in[1];
//     /// mix between dry and wet signal, clamped [0,1]
//     input& mix = ctrl.in[2];
//
//     array<input>& inputs() override { return in; }
//     array<output>& outputs() override { return out; }
//
//     void tick(T deltaTime) override
//     {
//       // delay time in samples, don't allow negative delay times because our buffer will fill up
//       T dts = clamp<T>(ctrl.in[0] / deltaTime, 1, MAX_BUFFER_SIZE);
//       T fbk = clamp<T>(ctrl.in[1], -1., 1.);
//       T fade = clamp<T>(ctrl.in[2], 0., 1.);
//       for (int i = 0; i < IO; ++i)
//       {
//         // audio input plus wet signal feedback from last tick
//         T dry = audio.in[i] + ctrl.out[i] * fbk;
//         if (buffer[i].capacity())
//         {
//           buffer[i].push(dry);
//         }
//
//         T wet = 0;
//
//         // NOTE: this method seems to work ok for fixed delay times and slowly changing delay times,
//         // but fast delay time changes generate gross artifacts,
//         // which will require a resampling solution to deal with (ala VCV Rack's Delay),
//         // meaning that this delay class is not entirely suitable for use in a flanger.
//         //
//         // how many samples do we need to pull out of the buffer to catch up with our delay time
//         T read = buffer[i].size() - dts;
//         if (read >= 0)
//         {
//           T s = buffer[i].front();
//           while (read >= 1)
//           {
//             s = buffer[i].pop();
//             read -= 1;
//           }
//
//           wet = s + (buffer[i].front() - s)*(1 - read);
//         }
//
//         ctrl.out[i] = wet;
//
//         audio.out[i] = dry + (wet - dry)*fade;
//       }
//     }
//
//   private:
//     ringbuffer<T, MAX_BUFFER_SIZE> buffer[IO];
//     io<T, IO, IO> audio;
//     io<T, 3,  IO> ctrl; // outputs here hold the wet signal from the previous tick
//   };
}

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
