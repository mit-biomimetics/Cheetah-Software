#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Utilities/BSplineBasic.h"
#include "Utilities/BezierCurve.h"
#include "Utilities/Utilities_print.h"
#include "Utilities/save_file.h"
#include "Utilities/utilities.h"
#include <Configuration.h>
#include "../user/WBC_Controller/WBC_States/Bounding/CtrlSet/ImpulseCurve.hpp"
#include "Utilities/Timer.h"

std::string folder_name = "/common/test/test_data/";

TEST(Spline, BezierCurve_test) {
  constexpr int dim = 3;
  constexpr int num_ctrl_pt = 4;

  BezierCurve<double, dim, num_ctrl_pt> bc;
  double** ctrl_pt = new double*[num_ctrl_pt];

  for (int i(0); i < 4; ++i) {
    ctrl_pt[i] = new double[dim];
  }
  double gap(0.1);
  // Initial
  ctrl_pt[0][0] = 0.0;
  ctrl_pt[0][1] = 0.0;
  ctrl_pt[0][2] = 0.4;

  // mid 1
  ctrl_pt[1][0] = 0.4;
  ctrl_pt[1][1] = -0.7;
  ctrl_pt[1][2] = ctrl_pt[0][2] + gap * 0.3;

  // mid 2
  ctrl_pt[2][0] = 2.5;
  ctrl_pt[2][1] = 0.5;
  ctrl_pt[2][2] = ctrl_pt[1][2] + gap * 0.6;

  // Final
  ctrl_pt[3][0] = 4.2;
  ctrl_pt[3][1] = -1.0;
  ctrl_pt[3][2] = ctrl_pt[2][2] + gap * 1.0;

  double end_time(3.);
  bc.SetParam(ctrl_pt, end_time);

  double curve_pt[dim];
  double curve_vel[dim];

  // 2nd curve
  double mid_time(2.);
  bc.getCurvePoint(mid_time, curve_pt);
  for (size_t i(0); i < dim; ++i) ctrl_pt[0][i] = curve_pt[i];

  ctrl_pt[1][2] = ctrl_pt[0][2] + gap;
  ctrl_pt[2][2] = ctrl_pt[1][2] + gap * 1.3;
  ctrl_pt[3][2] = ctrl_pt[2][2] + gap * 1.7;

  BezierCurve<double, dim, num_ctrl_pt> bc2;
  bc2.SetParam(ctrl_pt, end_time);

  // std::string folder_name = "/common/test/test_data/";
  //create_folder(folder_name);
  double t;
  double dt = 0.005;
  for (int i(0); i < 1001; ++i) {
    t = (double)i * dt;
    if (t < mid_time) {
      bc.getCurvePoint(t, curve_pt);
      bc.getCurveVelocity(t, curve_vel);
    } else {
      bc2.getCurvePoint(t - mid_time, curve_pt);
      bc2.getCurveVelocity(t - mid_time, curve_vel);
    }
//    saveVector(curve_pt, folder_name, "bezier_pos", dim);
//    saveVector(curve_vel, folder_name, "bezier_vel", dim);
//    saveValue(t, folder_name, "bz_time");
  }

  bc.getCurvePoint(0., curve_pt);
  EXPECT_TRUE(fpEqual(curve_pt[0], 0., .0001));
  EXPECT_TRUE(fpEqual(curve_pt[1], 0.0, .0001));
  EXPECT_TRUE(fpEqual(curve_pt[2], 0.4, .0001));

  bc.getCurvePoint(end_time, curve_pt);
  EXPECT_TRUE(fpEqual(curve_pt[0], 4.2, .0001));
  EXPECT_TRUE(fpEqual(curve_pt[1], -1.0, .0001));
  // EXPECT_TRUE(fpEqual(curve_pt[2], 0.9, .0001));

  for (int i(0); i < num_ctrl_pt; ++i) {
    delete[] ctrl_pt[i];
  }
  delete[] ctrl_pt;
}

TEST(Spline, BSpline_test) {
  constexpr int dim_bs = 3;
  constexpr int degree = 3;
  constexpr int num_middle_pt = 3;
  constexpr int ini_cstr_level = 2;
  constexpr int fin_cstr_level = 2;

  BS_Basic<double, dim_bs, degree, num_middle_pt, ini_cstr_level,
           fin_cstr_level>
      bs;
  double ini_pt[dim_bs * (degree)];
  double fin_pt[dim_bs * (degree)];

  for (int i(0); i < dim_bs * degree; ++i) {
    ini_pt[i] = 0.;
    fin_pt[i] = 0.;
  }
  double height_gap = 0.3;
  ini_pt[0] = 0.0;  // Initial
  ini_pt[1] = 0.0;
  ini_pt[2] = 0.4;

  double** middle_pt = NULL;
  if (num_middle_pt > 1) {
    middle_pt = new double*[num_middle_pt];
    for (int i(0); i < num_middle_pt; ++i) {
      middle_pt[i] = new double[dim_bs];
    }
    // Middle 1
    middle_pt[0][0] = 0.4;
    middle_pt[0][1] = -0.7;
    middle_pt[0][2] = ini_pt[2] + height_gap * 0.25;

    // Middle 2
    middle_pt[1][0] = 0.4;
    middle_pt[1][1] = -0.7;
    middle_pt[1][2] = ini_pt[2] + height_gap * 0.5;

    // Middle 3
    middle_pt[2][0] = 2.5;
    middle_pt[2][1] = 0.5;
    middle_pt[2][2] = ini_pt[2] + height_gap * 0.75;
  }
  // Final
  fin_pt[0] = 4.2;
  fin_pt[1] = -1.0;
  fin_pt[2] = ini_pt[2] + height_gap;

  fin_pt[3] = 0.0;  // Final vel
  fin_pt[4] = 0.0;
  fin_pt[5] = 0.0;

  fin_pt[6] = 0.0;  // Final acc
  fin_pt[7] = 0.0;
  fin_pt[8] = 0.0;

  double end_time(4.);
  bs.SetParam(ini_pt, fin_pt, middle_pt, end_time);

  double curve_pt[dim_bs];
  double curve_vel[dim_bs];
  double curve_acc[dim_bs];
  double curve_d3[dim_bs];

  // Spline 2
  double mid_time(1.);
  bs.getCurvePoint(mid_time, curve_pt);
  bs.getCurveDerPoint(mid_time, 1, curve_vel);
  bs.getCurveDerPoint(mid_time, 2, curve_acc);
  if (degree > 3) {
    bs.getCurveDerPoint(mid_time, 3, curve_d3);
  }
  BS_Basic<double, dim_bs, degree, num_middle_pt, ini_cstr_level,
           fin_cstr_level>
      bs_2;
  double ini_pt2[dim_bs * degree];
  double fin_pt2[dim_bs * degree];
  for (int i(0); i < dim_bs * degree; ++i) {
    ini_pt2[i] = 0.;
    fin_pt2[i] = 0.;
  }
  for (int i(0); i < 3; ++i) {
    ini_pt2[i] = curve_pt[i];  // Initial
    ini_pt2[i + 3] = curve_vel[i];
    ini_pt2[i + 6] = curve_acc[i];

    if (degree > 3) ini_pt2[i + 9] = curve_d3[i];
  }

  double** middle_pt2 = NULL;
  if (num_middle_pt > 0) {
    middle_pt2 = new double*[num_middle_pt];
    for (int i(0); i < num_middle_pt; ++i) {
      middle_pt2[i] = new double[dim_bs];
    }
    // Middle 1
    middle_pt2[0][0] = 0.4;
    middle_pt2[0][1] = -0.7;
    middle_pt2[0][2] = ini_pt2[2] + height_gap * 0.25;

    // Middle 2
    middle_pt2[1][0] = 0.4;
    middle_pt2[1][1] = -0.7;
    middle_pt2[1][2] = ini_pt2[2] + height_gap * 0.5;

    // Middle 3
    middle_pt2[2][0] = 2.5;
    middle_pt2[2][1] = 0.5;
    middle_pt2[2][2] = ini_pt2[2] + height_gap * 0.75;
  }

  // Final
  fin_pt2[0] = 4.2;
  fin_pt2[1] = -1.0;
  fin_pt2[2] = ini_pt2[2] + height_gap;

  bs_2.SetParam(ini_pt2, fin_pt2, middle_pt2, end_time);

  double t;
  double dt(0.005);
  for (int i(0); i < 1001; ++i) {
    t = (double)i * dt;

    if (t < mid_time) {
      bs.getCurvePoint(t, curve_pt);
      bs.getCurveDerPoint(t, 1, curve_vel);
      bs.getCurveDerPoint(t, 2, curve_acc);
    } else {
      bs_2.getCurvePoint(t - mid_time, curve_pt);
      bs_2.getCurveDerPoint(t - mid_time, 1, curve_vel);
      bs_2.getCurveDerPoint(t - mid_time, 2, curve_acc);
    }

    // saveVector(curve_pt, folder_name, "bspline_pos", dim_bs);
    // saveVector(curve_vel, folder_name, "bspline_vel", dim_bs);
    // saveVector(curve_acc, folder_name, "bspline_acc", dim_bs);
    // saveValue(t, folder_name, "bs_time");
  }

  bs.getCurvePoint(0., curve_pt);
  EXPECT_TRUE(fpEqual(curve_pt[0], 0., .0001));
  EXPECT_TRUE(fpEqual(curve_pt[1], 0.0, .0001));
  EXPECT_TRUE(fpEqual(curve_pt[2], 0.4, .0001));

  bs.getCurvePoint(end_time, curve_pt);
  EXPECT_TRUE(fpEqual(curve_pt[0], 4.2, .0001));
  EXPECT_TRUE(fpEqual(curve_pt[1], -1.0, .0001));
  // EXPECT_TRUE(fpEqual(curve_pt[2], ini_pt[2]+ height_gap, .0001));

  bs.getCurveDerPoint(end_time, 2, curve_acc);
  EXPECT_TRUE(fpEqual(curve_acc[0], 0., .0001));
  EXPECT_TRUE(fpEqual(curve_acc[1], 0.0, .0001));
  EXPECT_TRUE(fpEqual(curve_acc[2], 0.0, .0001));

  for (int i(0); i < num_middle_pt; ++i) {
    delete[] middle_pt[i];
  }
  delete[] middle_pt;
}

TEST(Spline, BSpline_1D_test) {
  constexpr int dim_bs = 1;
  constexpr int degree = 3;
  constexpr int num_middle_pt = 3;
  constexpr int ini_cstr_level = 2;
  constexpr int fin_cstr_level = 2;

  BS_Basic<double, dim_bs, degree, num_middle_pt, ini_cstr_level,
           fin_cstr_level>
      bs;
  double ini_pt[dim_bs * (degree)];
  double fin_pt[dim_bs * (degree)];

  for (int i(0); i < dim_bs * degree; ++i) {
    ini_pt[i] = 0.;
    fin_pt[i] = 0.;
  }
  double height_gap = 0.3;
  ini_pt[0] = 0.2;  // Initial
  ini_pt[1] = 0.0;
  ini_pt[2] = 0.0;

  double** middle_pt = NULL;
  if (num_middle_pt > 1) {
    middle_pt = new double*[num_middle_pt];
    for (int i(0); i < num_middle_pt; ++i) {
      middle_pt[i] = new double[dim_bs];
    }
    // Middle 1
    middle_pt[0][0] = ini_pt[0] + height_gap * 0.25;

    // Middle 2
    middle_pt[1][0] = ini_pt[0] + height_gap * 0.5;

    // Middle 3
    middle_pt[2][0] = ini_pt[0] + height_gap * 0.75;
  }
  // Final
  fin_pt[0] = ini_pt[0] + height_gap;

  double end_time(5.);
  bs.SetParam(ini_pt, fin_pt, middle_pt, end_time);

  double curve_pt[dim_bs];
  curve_pt[0] = 0.;
  double curve_vel[dim_bs];
  curve_vel[0] = 0.;
  double curve_acc[dim_bs];
  curve_acc[0] = 0.;

  double t;
  double dt(0.005);
  Vec3<double> bs_pva;
  for (int i(0); i < 1001; ++i) {
    t = (double)i * dt;

    bs.getCurvePoint(t, curve_pt);
    bs.getCurveDerPoint(t, 1, curve_vel);
    bs.getCurveDerPoint(t, 2, curve_acc);

    bs_pva[0] = curve_pt[0];
    bs_pva[1] = curve_vel[0];
    bs_pva[2] = curve_acc[0];
//    saveVector(bs_pva, folder_name, "bspline_1d");
//    saveValue(t, folder_name, "bs_time_1d");
  }

  for (int i(0); i < num_middle_pt; ++i) {
    delete[] middle_pt[i];
  }
  delete[] middle_pt;
}

TEST(Spline, ImpulseCurve) {
  ImpulseCurve<double> curve;

  double apex_value(2.9);
  double time(0.1);
  curve.setCurve(apex_value, time);

  double curve_pt;
  double t;
  double dt(0.001);
  double sum(0.);
  for (int i(-10); i < 501; ++i) {
    t = (double)i * dt;

    curve_pt = curve.getValue(t);
    sum += (curve_pt * dt);
//    saveValue(curve_pt, folder_name, "curve_pos");
//    saveValue(t, folder_name, "curve_time");
  }
  double integrated_value = 0.7 * apex_value * time;

  // printf("sum, integration: (%f, %f)\n", sum, integrated_value);

  EXPECT_TRUE(fpEqual(curve.getValue(time / 2.), apex_value, .0001));
  EXPECT_TRUE(fpEqual(integrated_value, sum, .0001));
}
