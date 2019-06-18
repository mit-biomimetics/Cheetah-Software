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

/// \cond INTERNAL

extern "C" {

  extern void sqic(
   const casadi_int *m, // Number of constraints + 1 (for the objective)
   const casadi_int* n, // Number of decision variables
   const casadi_int* nnzA, // Number of nonzeros in objective-augmented linear constraint matrix  A
   const casadi_int *indA, // colind of Compressed Column Storage A , length: nnzA
   const casadi_int *locA, // row of  Compressed Column Storage A, length n + 1
   const double *valA, // Values of A
   const double* bl, // Lower bounds to decision variables + objective
   const double* bu, // Upper bounds to decision variables + objective
   const casadi_int *hEtype, // ?
   const casadi_int *hs, // ?
   double *x,  // Decision variables + evaluated linear constraints ((initial + optimal), length n+m
   double *pi, // ?
   double *rc, // Multipliers (initial + optimal), length n+m
   const casadi_int* nnzH, // Number of nonzeros in full hessian H
   const casadi_int* indH, // colind of Compressed Column Storage H , length: nnzH
   const casadi_int* locH, // row of  Compressed Column Storage H, length n + 1
   double* valH);

  extern void sqicSolve(
   double* Obj); // Output: hessian part of the resulting objective

  extern void sqicSolveStabilized(
   double* Obj, // Output: hessian part of the resulting objective
   double *mu,
   casadi_int *lenpi,
   double* piE);

  extern void sqicDestroy();
}
/// \endcond
