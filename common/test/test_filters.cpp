/*! @file test_filters.cpp
 *  @brief Test filter functions
 *
 * Test the various filters
 */

#include "Dynamics/spatial.h"
#include "Math/FirstOrderIIRFilter.h"
#include "Math/orientation_tools.h"
#include "cppTypes.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace ori;
using namespace spatial;

/*!
 * Test the firstOrderIIR filter with Eigen Vec3's
 */
TEST(Filters, firstOrderIIR) {
  Vec3<double> x0(1, 2, 3);
  FirstOrderIIRFilter<Vec3<double>, double> filt(0.2, x0);
  Vec3<double> update(5, 5, 5);
  filt.update(update);
  Vec3<double> expected = x0 * 0.8 + 0.2 * update;
  EXPECT_TRUE(almostEqual(expected, filt.get(), .0001));
}
