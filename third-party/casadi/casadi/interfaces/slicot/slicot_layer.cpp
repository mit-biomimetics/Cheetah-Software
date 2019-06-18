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


#include "slicot_layer.hpp"

#include <vector>

// Need an 8-byte integer since libslicot0 is compiled with  fdefault-integer-8
typedef long long int f_int;

extern "C" {
  int mb03vd_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              double* tau, f_int* ldtau, double* dwork, f_int *info);
  int mb03vy_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              const double* tau, f_int* ldtau, double* dwork, f_int *ld_work, f_int *info);
  int mb03wd_(char* job, char* compz, f_int* n, f_int* p, f_int* ilo, f_int* ihi, f_int* iloz,
              f_int* ihiz, double *h, f_int* ldh1, f_int* ldh2, double* z, f_int* ldz1,
              f_int* ldz2, double* wr, double *wi, double* dwork, f_int *ld_work, f_int *info);


  int mb05nd_(f_int* n, double* delta, const double* a, f_int* lda,
                double* ex, f_int* ldex, double * exint, f_int* ldexin,
                double* tol, f_int* iwork, double * dwork, f_int* ldwork, f_int *info);

}

namespace casadi {

  int slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                     int ldtau, double * dwork) {
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ret_ = 0;

     mb03vd_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ret_);
     return ret_;
  }

  int slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                     const double * tau, int ldtau, double * dwork, int ldwork) {
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ldwork_ = ldwork;
     f_int ret_=0;
     mb03vy_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ldwork_, &ret_);

     return ret_;

  }

  int slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                     double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                     double *wi, double * dwork, int ldwork) {
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int iloz_ = ilo;
     f_int ihiz_ = ihi;
     f_int ldh1_ = ldh1;
     f_int ldh2_ = ldh2;
     f_int ldz1_ = ldz1;
     f_int ldz2_ = ldz2;
     f_int ldwork_ = ldwork;
     f_int ret_ = 0;
     mb03wd_(&job, &compz, &n_, &p_, &ilo_, &ihi_, &iloz_, &ihiz_, h, &ldh1_, &ldh2_,
             z, &ldz1_, &ldz2_, wr, wi, dwork, &ldwork_, &ret_);

     return ret_;
  }


  int slicot_mb05nd(int n, double delta, const double* a, int lda,
                     double* ex, int ldex, double * exint, int ldexin,
                     double tol, int* iwork, double * dwork, int ldwork) {
    f_int n_ = n;
    f_int lda_ = lda;
    f_int ldex_ = ldex;
    f_int ldexin_ = ldexin;
    f_int* iwork_ = reinterpret_cast<f_int*>(iwork);
    f_int ldwork_ = ldwork;
    f_int ret_ = 0;
    mb05nd_(&n_, &delta, a, &lda_, ex, &ldex_, exint, &ldexin_,
      &tol, iwork_, dwork, &ldwork_, &ret_);

    return ret_;
  }

} // namespace casadi
