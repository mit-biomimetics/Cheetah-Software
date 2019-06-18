/*! @file test_orientation_tools.cpp
 *  @brief Test orientation functions
 *
 * Test the orientation related functions in orientation_tools.h
 * Does not check any spatial stuff
 */


#include "cppTypes.h"
#include "Math/orientation_tools.h"
#include "Collision/collision_model.h"
#include "Dynamics/spatial.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

using namespace ori;
using namespace spatial;

/*!
 * Test the ground contact model
 */
TEST(Contact, normalContact) {
  Eigen::Matrix<double, 3, 10> pMATLAB, pdMATLAB, fMATLAB;
  Eigen::Matrix<double, 2, 10> uMATLAB, uOutMATLAB;

  pMATLAB << 0, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000,
          9.0000, 8.0000, 7.0000, 6.0000, 5.0000, 4.0000, 3.0000, 2.0000, 1.0000, 45.0000,
          0.1000, -0.1000, -0.0200, -0.0300, 0.0400, -1.0000, -2.0000, -0.0001, -0.1000, -0.2000;

  pdMATLAB <<
           1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000,
          -1.0000, -2.0000, -4.0000, -2.0000, 3.0000, 2.0000, 1.0000, 2.0000, 5.0000, 3.0000,
          -0.0010, -2.0000, -1.0000, 2.0000, -1.0000, -0.0010, 0.0020, 0.0200, -0.1000, 0.3000;

  uMATLAB << 0, 0, 0, 0, 0.2000, 0.2000, 0.3000, 0.4000, 0.5000, 0.2000,
          -0.2000, 0, 0, 0, -0.3000, 0.3000, 0.1000, 0.3000, 0.2000, 0;

  fMATLAB << 0, -15.8114, -10.6066, 0, 0, -68.0631, -201.0705, 0, -6.1473, 0,
          0, 15.8114, 14.1421, 0, 0, -23.9382, -29.4963, 0, -3.3903, 0,
          0, 16.1909, 3.5695, 0, 0, 12.0250, 33.8704, 0, 1.1700, 0;

  uOutMATLAB <<
             0, 2.0000, 3.0000, 0, 0.1040, 2.8265, 5.8431, 0.2080, 1.0376, 0.1040,
          -0.1040, -2.0000, -4.0000, 0, -0.1560, 1.1135, 0.8863, 0.1560, 0.5328, 0;


  vectorAligned<Vec3<double>> p, pd, f, fref;
  vectorAligned<Vec2<double>> u, uref;

  f.resize(10);

  for (size_t i = 0; i < 10; i++) {
    p.push_back(pMATLAB.block<3, 1>(0, i));
    pd.push_back(pdMATLAB.block<3, 1>(0, i));
    u.push_back(uMATLAB.block<2, 1>(0, i));
    fref.push_back(fMATLAB.block<3,1>(0,i));
    uref.push_back(uOutMATLAB.block<2,1>(0,i));
  }

  groundContactModel<double>(p, pd, u, f, 12, 25, 6, 1);

  for(size_t i = 0; i < 10; i++) {
    EXPECT_TRUE(almostEqual(fref[i], f[i], .0001));
    EXPECT_TRUE(almostEqual(uref[i], u[i], .0002));
  }
}

/*!
 * Test the ground contact model with offset set to identity
 */
TEST(Contact, normalContactWithNoOffset) {
  Eigen::Matrix<double, 3, 10> pMATLAB, pdMATLAB, fMATLAB;
  Eigen::Matrix<double, 2, 10> uMATLAB, uOutMATLAB;

  pMATLAB << 0, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000,
          9.0000, 8.0000, 7.0000, 6.0000, 5.0000, 4.0000, 3.0000, 2.0000, 1.0000, 45.0000,
          0.1000, -0.1000, -0.0200, -0.0300, 0.0400, -1.0000, -2.0000, -0.0001, -0.1000, -0.2000;

  pdMATLAB <<
           1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000,
          -1.0000, -2.0000, -4.0000, -2.0000, 3.0000, 2.0000, 1.0000, 2.0000, 5.0000, 3.0000,
          -0.0010, -2.0000, -1.0000, 2.0000, -1.0000, -0.0010, 0.0020, 0.0200, -0.1000, 0.3000;

  uMATLAB << 0, 0, 0, 0, 0.2000, 0.2000, 0.3000, 0.4000, 0.5000, 0.2000,
          -0.2000, 0, 0, 0, -0.3000, 0.3000, 0.1000, 0.3000, 0.2000, 0;

  fMATLAB << 0, -15.8114, -10.6066, 0, 0, -68.0631, -201.0705, 0, -6.1473, 0,
          0, 15.8114, 14.1421, 0, 0, -23.9382, -29.4963, 0, -3.3903, 0,
          0, 16.1909, 3.5695, 0, 0, 12.0250, 33.8704, 0, 1.1700, 0;

  uOutMATLAB <<
             0, 2.0000, 3.0000, 0, 0.1040, 2.8265, 5.8431, 0.2080, 1.0376, 0.1040,
          -0.1040, -2.0000, -4.0000, 0, -0.1560, 1.1135, 0.8863, 0.1560, 0.5328, 0;


  vectorAligned<Vec3<double>> p, pd, f, fref;
  vectorAligned<Vec2<double>> u, uref;

  f.resize(10);

  for (size_t i = 0; i < 10; i++) {
    p.push_back(pMATLAB.block<3, 1>(0, i));
    pd.push_back(pdMATLAB.block<3, 1>(0, i));
    u.push_back(uMATLAB.block<2, 1>(0, i));
    fref.push_back(fMATLAB.block<3,1>(0,i));
    uref.push_back(uOutMATLAB.block<2,1>(0,i));
  }

  Mat6<double> X = Mat6<double>::Identity();
  groundContactModelWithOffset<double>(p, pd, u, f, 12, 25, 6, 1, X);

  for(size_t i = 0; i < 10; i++) {
    EXPECT_TRUE(almostEqual(fref[i], f[i], .0001));
    EXPECT_TRUE(almostEqual(uref[i], u[i], .0002));
  }
}

/*!
 * Test the ground contact model with a negative z offset
 */
TEST(Contact, normalContactWithOffset) {
  Eigen::Matrix<double, 3, 10> pMATLAB, pdMATLAB, fMATLAB;
  Eigen::Matrix<double, 2, 10> uMATLAB, uOutMATLAB;

  pMATLAB << 0, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000,
          9.0000, 8.0000, 7.0000, 6.0000, 5.0000, 4.0000, 3.0000, 2.0000, 1.0000, 45.0000,
          0.1000 - 1, -0.1000 - 1, -0.0200 - 1, -0.0300 - 1, 0.0400 - 1, -1.0000 - 1, -2.0000 - 1, -0.0001 - 1, -0.1000 - 1, -0.2000 - 1;

  pdMATLAB <<
           1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000,
          -1.0000, -2.0000, -4.0000, -2.0000, 3.0000, 2.0000, 1.0000, 2.0000, 5.0000, 3.0000,
          -0.0010, -2.0000, -1.0000, 2.0000, -1.0000, -0.0010, 0.0020, 0.0200, -0.1000, 0.3000;

  uMATLAB << 0, 0, 0, 0, 0.2000, 0.2000, 0.3000, 0.4000, 0.5000, 0.2000,
          -0.2000, 0, 0, 0, -0.3000, 0.3000, 0.1000, 0.3000, 0.2000, 0;

  fMATLAB << 0, -15.8114, -10.6066, 0, 0, -68.0631, -201.0705, 0, -6.1473, 0,
          0, 15.8114, 14.1421, 0, 0, -23.9382, -29.4963, 0, -3.3903, 0,
          0, 16.1909, 3.5695, 0, 0, 12.0250, 33.8704, 0, 1.1700, 0;

  uOutMATLAB <<
             0, 2.0000, 3.0000, 0, 0.1040, 2.8265, 5.8431, 0.2080, 1.0376, 0.1040,
          -0.1040, -2.0000, -4.0000, 0, -0.1560, 1.1135, 0.8863, 0.1560, 0.5328, 0;


  vectorAligned<Vec3<double>> p, pd, f, fref;
  vectorAligned<Vec2<double>> u, uref;

  f.resize(10);

  for (size_t i = 0; i < 10; i++) {
    p.push_back(pMATLAB.block<3, 1>(0, i));
    pd.push_back(pdMATLAB.block<3, 1>(0, i));
    u.push_back(uMATLAB.block<2, 1>(0, i));
    fref.push_back(fMATLAB.block<3,1>(0,i));
    uref.push_back(uOutMATLAB.block<2,1>(0,i));
  }

  Mat6<double> X = createSXform(Mat3<double>::Identity(), Vec3<double>(0,0,-1));
  groundContactModelWithOffset<double>(p, pd, u, f, 12, 25, 6, 1, X);

  for(size_t i = 0; i < 10; i++) {
    EXPECT_TRUE(almostEqual(fref[i], f[i], .0001));
    EXPECT_TRUE(almostEqual(uref[i], u[i], .0002));
  }
}
