/*! @file test_math.cpp
 *  @brief Test math functions
 *
 */

#include "Math/Interpolation.h"
#include "Math/MathUtilities.h"
#include "cppTypes.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// test linear interpolation
TEST(Math, InterpolateLinear) {
  Vec3<double> y0(1, 2, 3);
  Vec3<double> y1(6, 5, 4);

  Vec3<double> yInterpRef = y0 + 0.75 * (y1 - y0);

  Vec3<double> yInterp = Interpolate::lerp(y0, y1, 0.75);

  EXPECT_TRUE(almostEqual(yInterp, yInterpRef, .001));
}

// test bezier interpolation derivatives
TEST(Math, BezierInterpolationFirstDerivative) {
  Vec3<double> y0(1, 2, 3);
  Vec3<double> y1(6, 5, 4);

  double dx = .00001;
  double x = 0.23;
  Vec3<double> yInterpX = Interpolate::cubicBezier(y0, y1, x);
  Vec3<double> yInterpXdx = Interpolate::cubicBezier(y0, y1, x + dx);

  Vec3<double> dyInterpDxRef = (yInterpXdx - yInterpX) / dx;

  Vec3<double> dyInterpDx = Interpolate::cubicBezierFirstDerivative(y0, y1, x);

  EXPECT_TRUE(almostEqual(dyInterpDx, dyInterpDxRef, .001));

  for (int i = 0; i < 3; i++) {
    EXPECT_TRUE(yInterpX(i) > y0(i) && yInterpX(i) < y1(i));
  }
}

// test bezier interpolation derivatives
TEST(Math, BezierInterpolationSecondDerivative) {
  Vec3<double> y0(1, 2, 3);
  Vec3<double> y1(6, 5, 4);

  double dx = .00001;
  double x = 0.23;
  Vec3<double> yInterpX = Interpolate::cubicBezierFirstDerivative(y0, y1, x);
  Vec3<double> yInterpXdx =
      Interpolate::cubicBezierFirstDerivative(y0, y1, x + dx);

  Vec3<double> dyInterpDxRef = (yInterpXdx - yInterpX) / dx;

  Vec3<double> dyInterpDx = Interpolate::cubicBezierSecondDerivative(y0, y1, x);

  EXPECT_TRUE(almostEqual(dyInterpDx, dyInterpDxRef, .001));
}