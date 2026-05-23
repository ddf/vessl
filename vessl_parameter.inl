#pragma once

namespace vessl
{
template <size_t N>
VESSL_INLINE constexpr parameter::desc parameter::desc_list<N>::operator[](id_t id) const
{
  for (size_t i = 0; i < size; ++i)
  {
    if (list[i].id == id)
    {
      return list[i];
    }
  }
  return desc::empty();
}

template <typename T>
VESSL_INLINE T parameter::read() const
{
  switch (desc_.type)
  {
  case value_type::none: return T(0LL);
  case value_type::binary: return cast<T>(*static_cast<binary_t*>(data_));
  case value_type::digital: return cast<T>(*static_cast<digital_t*>(data_));
  case value_type::analog: return cast<T>(*static_cast<analog_t*>(data_));
  case value_type::phase: return cast<T>(*static_cast<phase_t*>(data_));
  case value_type::user: break;
  }
  // attempt to cast data to T, might work?
  return *static_cast<T*>(data_);
}

template <typename T>
VESSL_INLINE parameter & parameter::write(const T &value)
{
  switch (desc_.type)
  {
  case value_type::none: break;
  case value_type::binary: *static_cast<binary_t*>(data_)   = cast<binary_t>(value); break;
  case value_type::digital: *static_cast<digital_t*>(data_) = cast<digital_t>(value); break;
  case value_type::analog: *static_cast<analog_t*>(data_)   = cast<analog_t>(value); break;
  case value_type::phase: *static_cast<phase_t*>(data_)     = cast<phase_t>(value); break;
  // attempt to cast data to T, might work?
  case value_type::user: *static_cast<T*>(data_) = value; break;
  }
  return *this;
}

VESSL_INLINE parameter & parameter::operator=(const parameter &rhs)
{
  if (&rhs != this)
  {
    switch (desc_.type)
    {
    case value_type::none: break;
    case value_type::binary: *static_cast<binary_t*>(data_)   = rhs.read_binary(); break;
    case value_type::digital: *static_cast<digital_t*>(data_) = rhs.read_digital(); break;
    case value_type::analog: *static_cast<analog_t*>(data_)   = rhs.read_analog(); break;
    case value_type::phase: *static_cast<phase_t*>(data_)     = rhs.read_phase(); break;
    case value_type::user: 
      switch (rhs.desc_.type)
      {
    case value_type::none: break;
    case value_type::binary: write(rhs.read_binary()); break;
    case value_type::digital: write(rhs.read_digital()); break;
    case value_type::analog: write(rhs.read_analog()); break;
    case value_type::phase: write(rhs.read_phase()); break;
    case value_type::user: VASSERT(false, "Can't assign user parameter to user parameter with operator="); break;
      }
      break;
    }
  }
  return *this;
}

template<>
VESSL_INLINE constexpr parameter::value_type parameter::type_of<void*>() { return value_type::none; }

template<>
VESSL_INLINE constexpr parameter::value_type parameter::type_of<analog_t>() { return value_type::analog; }

template<>
VESSL_INLINE constexpr parameter::value_type parameter::type_of<digital_t>() { return value_type::digital; }

template<>
VESSL_INLINE constexpr parameter::value_type parameter::type_of<phase_t>() { return value_type::phase; }
  
template<>
struct parameter::data<void*>
{
  void* value = nullptr;
};
  
template<>
struct parameter::data<analog_t>
{
  analog_t value = 0.f;
};
  
template<>
struct parameter::data<digital_t>
{
  digital_t value = 0;
};
  
template<>
struct parameter::data<binary_t>
{
  binary_t value = false;
};
  
template<>
struct parameter::data<phase_t>
{
  phase_t value = phase_zero;
};

template<>
struct parameter::data<gain_t>
{
  gain_t value = gain_t();
};

template<>
struct parameter::data<duration_t>
{
  duration_t value = duration_t();
};

template <typename T>
VESSL_INLINE constexpr parameter param<T>::operator()(const char_t *name, parameter::id_t id) const
{
  return parameter(parameter::desc(name, id, parameter::type_of<T>()), *this);
}

VESSL_INLINE parameter parameter::none()
{
  data<void*> v; return parameter(desc::empty(), v);
}

} // namespace vessl