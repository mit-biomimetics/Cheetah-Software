/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_SLICOT_LA_HPP
#define CASADI_SLICOT_LA_HPP

/// \cond INTERNAL
namespace casadi {


  inline void dense_kron_stride(casadi_int n, casadi_int m,
      const double *A, const double *B, double *C,
      casadi_int strideA, casadi_int strideB, casadi_int strideC) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<n;++j) {
        C[strideC*i + j] = -A[strideA*(i/m) + (j/m)]*B[strideB*(i%m) + (j%m)];
      }
    }
  }

  inline void dense_mul_nt_stride(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C,
      casadi_int strideA, casadi_int strideB, casadi_int strideC) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        for (casadi_int k=0;k<l;++k) {
          C[strideC*i + j] += A[strideA*i + k]*B[strideB*j + k];
        }
      }
    }
  }

  //  A : n-by-l   B: m-by-l
  //  C = A*B'
  inline void dense_mul_nt(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C) {
    dense_mul_nt_stride(n, m, l, A, B, C, n, m, n);
  }

  inline void dense_mul_nn_stride(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C,
      casadi_int strideA, casadi_int strideB, casadi_int strideC) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        for (casadi_int k=0;k<l;++k) {
          C[strideC*i + j] += A[strideA*i + k]*B[strideB*k + j];
        }
      }
    }
  }

  inline void dense_copy_stride(casadi_int n, casadi_int m, const double *A, double *B,
      casadi_int strideA, casadi_int strideB) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        B[strideB*i + j] = A[strideA*i+j];
      }
    }
  }

  inline void dense_copy_t_stride(casadi_int n, casadi_int m, const double *A, double *B,
      casadi_int strideA, casadi_int strideB) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        B[strideB*j + i] = A[strideA*i+j];
      }
    }
  }

  //  A : n-by-l   B: l-by-m
  //  C = A*B
  inline void dense_mul_nn(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C) {
    dense_mul_nn_stride(n, m, l, A, B, C, n, l, n);
  }

  //  A : n-by-l   B: l-by-m
  //  C = A*B
  inline void dense_mul_nn2(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        for (casadi_int k=0;k<l;++k) {
          C[i + n*j] += A[i + n*k]*B[k + l*j];
        }
      }
    }
    //dense_mul_nn_stride(m, n, l, B, A, C, l, n, n);
  }

  //  A : l-by-n   B: l-by-m
  //  C = A'*B
  inline void dense_mul_tn(casadi_int n, casadi_int m, casadi_int l,
      const double *A, const double *B, double *C) {
    for (casadi_int i=0;i<n;++i) {
      for (casadi_int j=0;j<m;++j) {
        for (casadi_int k=0;k<l;++k) {
          C[n*i + j] += A[l*k + i]*B[l*k + j];
        }
      }
    }
  }

} // namespace casadi

/// \endcond
#endif // CASADI_SLICOT_LA_HPP
