////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2026 Damien Quartz
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

// Template specializations that use the optimized CMSIS library:
// http://www.keil.com/pack/doc/CMSIS/General/html/index.html
// #define ARM_CORTEX before including vessl.h to enable these.
#include <arm_math.h>

namespace vessl
{
  template<>
  inline void array<float32_t>::copyTo(array dest)
  {
    VASSERT(this->getSize() <= dest.getSize(), "Not enough room in destination for this array");
    arm_copy_f32(data, dest.getData(), this->getSize());
  }

  template<>
  inline void array<float32_t>::fill(float value)
  {
    arm_fill_f32(value, data, size);  
  }
  
  template<>
  inline array<float32_t> array<float32_t>::offset(float value, array dest) const
  {
    VASSERT(this->getSize() <= dest.getSize(), "Not enough room in destination for this array");
    arm_offset_f32(data, value, dest.getData(), this->getSize());
    return dest;
  }
  
  template<>
  inline array<float32_t> array<float32_t>::add(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_add_f32(data, other.data, dest.data, size);
    return dest;
  }

  template<>
  inline array<float32_t> array<float32_t>::subtract(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_sub_f32(data, other.data, dest.data, size);
    return dest;
  }

  template<>
  inline array<float32_t> array<float32_t>::scale(float32_t value, array dest) const
  {
    VASSERT(size <= dest.size, "Destination array is not large enough");
    arm_scale_f32(data, value, dest.data, size);
    return dest;
  }

  template<>
  inline array<float32_t> array<float32_t>::multiply(array other, array dest) const
  {
    VASSERT(size == other.size && size <= dest.size, "Arrays are different sizes or dest is not large enough.");
    arm_mult_f32(data, other.data, dest.data, size);
    return dest;
  }
  
  template<>
  struct matrixData<float32_t>
  {
    arm_matrix_instance_f32 inst;
    matrixData() { arm_mat_init_f32(&inst, 0, 0, nullptr); };
    matrixData(float32_t* data, uint32_t r, uint32_t c) { arm_mat_init_f32(&inst, r, c, data); }
    
    float32_t* operator*() { return inst.pData; }
    float32_t* operator*() const { return inst.pData; }
    [[nodiscard]] uint32_t rows() const { return inst.numRows; }
    [[nodiscard]] uint32_t cols() const { return inst.numCols; }
  };
  
  template<>
  inline matrix<float32_t> matrix<float32_t>::add(matrix other, matrix dest) const
  {
    VASSERT(rows == other.rows && rows == dest.rows && columns == other.columns && colums == dest.columns, "matrices do not have the same dimentions");
    arm_mat_add_f32(&data.inst, &other.data.inst, &dest.data.inst);
    return dest;
  }

  template<>
  inline matrix<float32_t> matrix<float32_t>::subtract(matrix other, matrix dest) const
  {
    VASSERT(rows == other.rows && rows == dest.rows && columns == other.columns && colums == dest.columns, "matrices do not have the same dimentions");
    arm_mat_sub_f32(&data.inst, &other.data.inst, &dest.data.inst);
    return dest;
  }
  
  template<>
  inline matrix<float32_t> matrix<float32_t>::scale(float32_t value, matrix dest) const
  {
    arm_mat_scale_f32(&data.inst, value, &dest.data.inst);
    return dest;
  }
  
  template<>
  inline matrix<float32_t> matrix<float32_t>::multiply(matrix other, matrix dest) const
  {
    VASSERT(getColumns() == other.getRows(), "Incompatible matrix sizes in operands");
    VASSERT(dest.getRows() == getRows(), "Incorrect number of rows in destination");
    VASSERT(dest.getColumns() == other.getColumns(), "Incorrect number of columns in destination");
    arm_mat_mult_f32(&data.inst, &other.data.inst, &dest.data.inst);
    return dest;
  }
  
  // @todo can implement when/if I update the version of CMSIS used by OWL
  // template<>
  // inline array<float32_t> matrix<float32_t>::multiply(array<float32_t> vector, array<float32_t> dest) const
  // {
  //   VASSERT(getColumns() == vector.getSize(), "Incompatible operands");
  //   VASSERT(dest.getSize() == getRows(), "Incompatible destination size");
  //   return dest;
  // }
  
  namespace math
  {
    template<>
    inline float32_t sin<float32_t, float32_t>(float32_t r) { return arm_sin_f32(r); }

    template<>
    inline float32_t sqrt(float32_t x)
    {
      float out;
      if (ARM_MATH_SUCCESS == arm_sqrt_f32(x, &out))
      {
        return out;
      }
      return 0;
    }
  }
  
  namespace filtering
  {
    template<size_t STAGES>
    template<class CoGen>
    struct biquad<STAGES>::df2T<float, CoGen> final : cascade<float, 2>
    {
      arm_biquad_cascade_df2T_instance_f32 inst;
          
      using biquad<STAGES>::cascade<float, 2>::coeff;
      using biquad<STAGES>::cascade<float, 2>::state;
      using biquad<STAGES>::cascade<float, 2>::getCoeffSize;
      using biquad<STAGES>::cascade<float, 2>::getStateSize;
    
      static CoGen cg;
  
      df2T() { arm_biquad_cascade_df2T_init_f32(&inst, STAGES, coeff.getData(), state.getData()); }
    
      size_t getStageCount() const { return STAGES; }
  
      void process(const float* source, float* dest, size_t blockSize, const args& args)
      {
        cg(coeff.getData(), args.omega(), args.q, args.g);
        arm_biquad_cascade_df2T_f32(&inst, source, dest, blockSize);
      }
    };
  }
}
