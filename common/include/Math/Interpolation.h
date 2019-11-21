/*! @file Interpolation.h
 *  @brief Utility functions to interpolate between two values
 *
 */

#ifndef PROJECT_INTERPOLATION_H
#define PROJECT_INTERPOLATION_H

#include <assert.h>
#include <type_traits>

namespace Interpolate {

/*!
 * Linear interpolation between y0 and yf.  x is between 0 and 1
 */
template <typename y_t, typename x_t>
y_t lerp(y_t y0, y_t yf, x_t x) {
  static_assert(std::is_floating_point<x_t>::value,
                "must use floating point value");
  assert(x >= 0 && x <= 1);
  return y0 + (yf - y0) * x;
}

/*!
 * Cubic bezier interpolation between y0 and yf.  x is between 0 and 1
 */
template <typename y_t, typename x_t>
y_t cubicBezier(y_t y0, y_t yf, x_t x) {
  static_assert(std::is_floating_point<x_t>::value,
                "must use floating point value");
  assert(x >= 0 && x <= 1);
  y_t yDiff = yf - y0;
  x_t bezier = x * x * x + x_t(3) * (x * x * (x_t(1) - x));
  return y0 + bezier * yDiff;
}

/*!
 * Cubic bezier interpolation derivative between y0 and yf.  x is between 0 and
 * 1
 */
template <typename y_t, typename x_t>
y_t cubicBezierFirstDerivative(y_t y0, y_t yf, x_t x) {
  static_assert(std::is_floating_point<x_t>::value,
                "must use floating point value");
  assert(x >= 0 && x <= 1);
  y_t yDiff = yf - y0;
  x_t bezier = x_t(6) * x * (x_t(1) - x);
  return bezier * yDiff;
}

/*!
 * Cubic bezier interpolation derivative between y0 and yf.  x is between 0 and
 * 1
 */
template <typename y_t, typename x_t>
y_t cubicBezierSecondDerivative(y_t y0, y_t yf, x_t x) {
  static_assert(std::is_floating_point<x_t>::value,
                "must use floating point value");
  assert(x >= 0 && x <= 1);
  y_t yDiff = yf - y0;
  x_t bezier = x_t(6) - x_t(12) * x;
  return bezier * yDiff;
}

}  // namespace Interpolate

#endif  // PROJECT_INTERPOLATION_H
