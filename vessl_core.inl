#pragma once

namespace vessl
{

// phase t <--> digital_t (@todo convert to degrees?)
template<>
VESSL_INLINE constexpr digital_t cast<digital_t, phase_t>(const phase_t& from)
{
  return from > phase_180 ? 1 : 0;
}

template<>
VESSL_INLINE constexpr phase_t cast<phase_t, digital_t>(const digital_t& from)
{
  return from <= 0 ? phase_zero : phase_360;
}

template<>
VESSL_INLINE constexpr analog_t cast<analog_t, phase_t>(const phase_t& from)
{
  return static_cast<analog_t>(from) / 4294967295.0f;
}

template<>
VESSL_INLINE constexpr phase_t cast<phase_t, analog_t>(const analog_t& from)
{
  return static_cast<phase_t>(from * 4294967295.0);
}

template<>
VESSL_INLINE constexpr phase_t cast<phase_t, size_t>(const size_t& from)
{
  return from % phase_360;
}
  
// phase_t <--> binary_t
template<>
VESSL_INLINE constexpr binary_t cast<binary_t, phase_t>(const phase_t& from)
{
  return from > phase_180;
}
  
template<>
VESSL_INLINE constexpr phase_t cast<phase_t, binary_t>(const binary_t& from)
{
  return from ? phase_360 : phase_zero;
}


template<typename T>
VESSL_INLINE matrix<T> matrix<T>::add(matrix other, matrix dest) const
{
  VASSERT(rows() == other.rows() 
    && rows() == dest.rows() 
    && columns() == other.columns() 
    && columns() == dest.columns()
    , "matrices do not have the same dimensions");
  
  array<T> lhs(data(), size());
  array<T> rhs(other.data(), other.size());
  array<T> dst(dest.data(), dest.size());
  lhs.add(rhs, dst);
  return dest;
}

template<typename T>
VESSL_INLINE matrix<T> matrix<T>::subtract(matrix other, matrix dest) const
{
  VASSERT(rows() == other.rows() 
    && rows() == dest.rows() 
    && columns() == other.columns() 
    && columns() == dest.columns()
    , "matrices do not have the same dimensions");
  
  array<T> lhs(data(), size());
  array<T> rhs(other.data(), other.size());
  array<T> dst(dest.data(), dest.size());
  lhs.subtract(rhs, dst);
  return dest;
}

template<typename T>
VESSL_INLINE matrix<T> matrix<T>::scale(T value, matrix dest) const
{
  array<T> lhs(data(), size());
  array<T> dst(dest.data(), dest.size());
  lhs.scale(value, dst);
  return dest;
}

template<typename T>
VESSL_INLINE matrix<T> matrix<T>::multiply(matrix other, matrix dest) const
{
  VASSERT(columns() == other.rows(), "Incompatible matrix sizes in operands");
  VASSERT(dest.rows() == rows(), "Incorrect number of rows in destination");
  VASSERT(dest.columns() == other.columns(), "Incorrect number of columns in destination");
  
  for(size_t i = 0; i < rows(); i++)
  {
    for(size_t j = 0; j < other.columns(); j++)
    {
      T accum = T(0LL);
      for(size_t k = 0; k < other.rows(); k++)
      {
        accum += get(i, k) * other.get(k, j);
      }
      dest.set(i, j, accum);
    }
  }
  return dest;
}

template<typename T>
VESSL_INLINE array<T> matrix<T>::multiply(const array<T>& vector, array<T> dest) const
{
  VASSERT(columns() == vector.size(), "Incompatible operands");
  VASSERT(dest.size() == rows(), "Incompatible destination size");
  
  for(size_t i = 0; i < rows(); i++)
  {
    T accum = 0;
    for(size_t j = 0; j < vector.size(); j++)
    {
      accum += get(i, j) * vector[j];
    }
    dest[i] = accum;
  }
  return dest;
}

template<typename I, typename O>
VESSL_INLINE void processor<I,O>::process(source<I>& in, sink<O>& out)
{
  while (in && out)
  {
    out.write(process(in.read()));
  }
}

template<typename I, typename O>
VESSL_INLINE void processor<I,O>::process(array<I> in, array<O> out)
{
  auto r = in.make_reader();
  auto w = out.make_writer();
  process(r, w);
}


template <typename T>
VESSL_INLINE void matrix<T>::clear()
{
  array<T> arr(data(), size()); 
  arr.fill(T(0LL));
}

// @todo ARM specialization
template<typename T>
VESSL_INLINE void array<T>::writer::write(const reader& r)
{
  size_t rsz = r.available();
  VASSERT(available() >= rsz, "Not enough space in writer for the contents of reader");
  const T* rh = *r;
  memcpy(static_cast<void*>(head_), static_cast<const void*>(rh), rsz * sizeof(T));
  head_ += rsz;
}
  
template<typename T>
VESSL_INLINE void array<T>::copy_to(array dest) const
{
  writer w(dest);
  reader r(*this);
  w.write(r);
}

template<typename T>
VESSL_INLINE void array<T>::fill(T value)
{
  writer w (*this);
  while (w)
  {
    w << value;
  }
}

template<typename T>
VESSL_INLINE array<T> array<T>::offset(T value, array dest) const
{
  VASSERT(size_ <= dest.size, "arrays are have different lengths or destination is too small");
  reader a(data_, size_);
  writer b(dest.data, dest.size);
  while (a)
  {
    b << a.read() + value;
  }
  return dest;
}

template<typename T>
VESSL_INLINE array<T> array<T>::add(array other, array dest) const
{
  VASSERT(size_ == other.size_
    && size_ <= dest.size_
    , "arrays are have different lengths or destination is too small");
  reader a(data_, size_);
  reader b(other.data_, other.size_);
  writer c(dest.data_, dest.size_);
  while (a)
  {
    c << a.read() + b.read();
  }
  return dest;
}

template<typename T>
VESSL_INLINE array<T> array<T>::subtract(array other, array dest) const
{
  VASSERT(size_ == other.size_ 
    && size_ <= dest.size_
    , "arrays are have different lengths or destination is too small");
  reader a(data_, size_);
  reader b(other.data_, other.size_);
  writer c(dest.data_, dest.size_);
  while (a)
  {
    c << a.read() - b.read();
  }
  return dest;
}

template<typename T>
VESSL_INLINE array<T> array<T>::scale(T value, array dest) const
{
  VASSERT(size_ <= dest.size_, "destination size is too small");
  reader a(data_, size_);
  writer b(dest.data_, dest.size_);
  while (a)
  {
    b << a.read() * value;
  }
  return dest;
}

template<typename T>
VESSL_INLINE array<T> array<T>::multiply(array other, array dest) const
{
  VASSERT(size_ == other.size_
    && size_ <= dest.size_
    , "arrays are have different lengths or destination is too small");
  reader a(data_, size_);
  reader b(other.data_, other.size_);
  writer c(dest.data_, dest.size_);
  while (a)
  {
    c << a.read() * b.read();
  }
  return dest;
}

template<typename T>
VESSL_INLINE typename array<T>::writer& operator<<(typename array<T>::writer& w, const T& v)
{
  w.write(v);
  return w;
}
} // namespace vessl