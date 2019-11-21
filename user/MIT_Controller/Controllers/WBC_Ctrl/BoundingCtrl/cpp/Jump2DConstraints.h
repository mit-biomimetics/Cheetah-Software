//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Jump2DConstraints.h
//
// MATLAB Coder version            : 3.2
// C/C++ source code generated on  : 26-Aug-2019 18:20:14
//
#ifndef JUMP2DCONSTRAINTS_H
#define JUMP2DCONSTRAINTS_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "Jump2DBounds_types.h"

// Function Declarations
extern void Jump2DConstraints(const double in1[6], const double in2[4], const
  double in3[6], const double in4[2], double dt, const double in6[4], double m,
  double Iyy, double g, double mu_g, double constraints[10]);

#endif

//
// File trailer for Jump2DConstraints.h
//
// [EOF]
//
