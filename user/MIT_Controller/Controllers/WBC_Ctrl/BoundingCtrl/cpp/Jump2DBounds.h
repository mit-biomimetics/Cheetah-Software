//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Jump2DBounds.h
//
// MATLAB Coder version            : 3.2
// C/C++ source code generated on  : 26-Aug-2019 18:20:14
//
#ifndef JUMP2DBOUNDS_H
#define JUMP2DBOUNDS_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "Jump2DBounds_types.h"

// Function Declarations
extern void Jump2DBounds(const double in1[6], const double in2[6], const double
  in3[4], const double in4[4], const double in5[2], double decision_vars_lb[10],
  double decision_vars_ub[10], double constraints_ub[10], double constraints_lb
  [10]);

#endif

//
// File trailer for Jump2DBounds.h
//
// [EOF]
//
