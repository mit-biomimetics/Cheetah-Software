/*! @file test_ControlParameters.cpp
 *  @brief
 *
 */

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "ControlParameters/ControlParameters.h"
#include "ControlParameters/RobotParameters.h"
#include "ControlParameters/SimulatorParameters.h"
#include "Math/MathUtilities.h"

class TestControlParameters : public ControlParameters {
 public:
  TestControlParameters()
      : ControlParameters("test-parameters"),
        test_double_param("test_double", test_double, collection),
        test_float_param("test_float", test_float, collection),
        test_integer_param("test_integer", test_integer, collection) {}

  double test_double;
  ControlParameter test_double_param;

  float test_float;
  ControlParameter test_float_param;

  s64 test_integer;
  ControlParameter test_integer_param;
};

class TestVectorControlParameters : public ControlParameters {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  TestVectorControlParameters()
      : ControlParameters("vector-params"),
        test_double_param("test_double", test_double, collection),
        test_float_param("test_float", test_float, collection) {}

  Vec3<double> test_double;
  ControlParameter test_double_param;

  Vec3<float> test_float;
  ControlParameter test_float_param;
};

TEST(ControlParams, testSimple) {
  TestControlParameters settings;

  EXPECT_FALSE(settings.isFullyInitialized());

  EXPECT_THROW(settings.initializeInteger("test_double", 1),
               std::runtime_error);
  EXPECT_THROW(settings.initializeDouble("not-a-real-thing", 1.2),
               std::runtime_error);

  settings.initializeDouble("test_double", 1.2);

  EXPECT_FALSE(settings.isFullyInitialized());
  EXPECT_TRUE(1.2 == *settings.collection.lookup("test_double")._value.d);
  EXPECT_TRUE(ControlParameterValueKind::DOUBLE ==
              settings.collection.lookup("test_double")._kind);
  EXPECT_TRUE("test_double" == settings.collection.lookup("test_double")._name);
  EXPECT_TRUE(settings.collection.lookup("test_double")._set);

  settings.initializeFloat("test_float", 1);
  EXPECT_FALSE(settings.isFullyInitialized());
  EXPECT_TRUE(1 == *settings.collection.lookup("test_float")._value.f);
  EXPECT_TRUE(ControlParameterValueKind::FLOAT ==
              settings.collection.lookup("test_float")._kind);
  EXPECT_TRUE("test_float" == settings.collection.lookup("test_float")._name);
  EXPECT_TRUE(settings.collection.lookup("test_float")._set);

  settings.initializeInteger("test_integer", 8);
  EXPECT_TRUE(settings.isFullyInitialized());
  EXPECT_TRUE(8 == *settings.collection.lookup("test_integer")._value.i);
  EXPECT_TRUE(ControlParameterValueKind::S64 ==
              settings.collection.lookup("test_integer")._kind);
  EXPECT_TRUE("test_integer" ==
              settings.collection.lookup("test_integer")._name);
  EXPECT_TRUE(settings.collection.lookup("test_integer")._set);
}

TEST(ControlParams, testIni) {
  // create and initialize some parameters
  TestControlParameters settings;
  EXPECT_FALSE(settings.isFullyInitialized());
  settings.initializeDouble("test_double", 2e-9);
  settings.initializeFloat("test_float", 1);
  settings.initializeInteger("test_integer", 8);
  EXPECT_TRUE(settings.isFullyInitialized());

  // write to ini:
  settings.writeToIniFile("control-params-test-file.ini");

  // read from ini
  TestControlParameters settingsFromIni;
  EXPECT_FALSE(settingsFromIni.isFullyInitialized());
  settingsFromIni.initializeFromIniFile("control-params-test-file.ini");
  printf("%s\n", settingsFromIni.generateUnitializedList().c_str());
  EXPECT_TRUE(settingsFromIni.isFullyInitialized());

  EXPECT_TRUE(
      fpEqual(settings.test_double, settingsFromIni.test_double, 1e-10));
  EXPECT_TRUE(fpEqual(settings.test_float, settingsFromIni.test_float, 1e-5f));
  EXPECT_TRUE(settings.test_integer == settingsFromIni.test_integer);
}

TEST(ControlParams, testVector) {
  TestVectorControlParameters settings;
  EXPECT_FALSE(settings.isFullyInitialized());
  Vec3<float> v3f(1, 2.2, 54.45);
  Vec3<double> v3d(.01, 1e-3, -2.23e-7);
  settings.initializeVec3f("test_float", v3f);
  settings.initializeVec3d("test_double", v3d);
  EXPECT_TRUE(settings.isFullyInitialized());

  settings.writeToYamlFile("test-yaml.yaml");

  TestVectorControlParameters settingsFromFile;
  EXPECT_FALSE(settingsFromFile.isFullyInitialized());
  settingsFromFile.initializeFromYamlFile("test-yaml.yaml");
  EXPECT_TRUE(settingsFromFile.isFullyInitialized());

  EXPECT_TRUE(almostEqual(v3f, settingsFromFile.test_float, 1e-6f));
  EXPECT_TRUE(almostEqual(v3d, settingsFromFile.test_double, 1e-10));
}

TEST(ControlParams, testYaml) {
  // create and initialize some parameters
  TestControlParameters settings;
  EXPECT_FALSE(settings.isFullyInitialized());
  settings.initializeDouble("test_double", 2e-9);
  settings.initializeFloat("test_float", 1);
  settings.initializeInteger("test_integer", 8);
  EXPECT_TRUE(settings.isFullyInitialized());

  // write to yaml
  settings.writeToYamlFile("control-params-test-file.yaml");

  // read from yaml
  TestControlParameters settingsFromYaml;
  EXPECT_FALSE(settingsFromYaml.isFullyInitialized());
  settingsFromYaml.initializeFromYamlFile("control-params-test-file.yaml");
  EXPECT_TRUE(settingsFromYaml.isFullyInitialized());

  EXPECT_TRUE(
      fpEqual(settings.test_double, settingsFromYaml.test_double, 1e-10));
  EXPECT_TRUE(fpEqual(settings.test_float, settingsFromYaml.test_float, 1e-5f));
  EXPECT_TRUE(settings.test_integer == settingsFromYaml.test_integer);
}

// check to see that the simulator default settings file contains all the
// simulator settings.
TEST(ControlParams, CheckSimulatorDefaults) {
  SimulatorControlParameters simParams;
  simParams.initializeFromYamlFile(getConfigDirectoryPath() +
                                   SIMULATOR_DEFAULT_PARAMETERS);
  if (!simParams.isFullyInitialized()) {
    printf("Missing parameters:\n%s\n",
           simParams.generateUnitializedList().c_str());
  }
  EXPECT_TRUE(simParams.isFullyInitialized());
}

// check to see that the simulator default settings file contains all the
// simulator settings.
TEST(ControlParams, CheckMiniCheetahDefaults) {
  RobotControlParameters robotParams;
  robotParams.initializeFromYamlFile(getConfigDirectoryPath() +
                                     MINI_CHEETAH_DEFAULT_PARAMETERS);
  if (!robotParams.isFullyInitialized()) {
    printf("Missing parameters:\n%s\n",
           robotParams.generateUnitializedList().c_str());
  }
  EXPECT_TRUE(robotParams.isFullyInitialized());
}

// check to see that the simulator default settings file contains all the
// simulator settings.
TEST(ControlParams, CheckCheetah3Defulats) {
  RobotControlParameters robotParams;
  robotParams.initializeFromYamlFile(getConfigDirectoryPath() +
                                     CHEETAH_3_DEFAULT_PARAMETERS);
  if (!robotParams.isFullyInitialized()) {
    printf("Missing parameters:\n%s\n",
           robotParams.generateUnitializedList().c_str());
  }
  EXPECT_TRUE(robotParams.isFullyInitialized());
}


TEST(ControlParams, CheckTypeRecognition) {
  EXPECT_TRUE(getControlParameterValueKindFromString("3.f") == ControlParameterValueKind::FLOAT);
  EXPECT_TRUE(getControlParameterValueKindFromString("3.") == ControlParameterValueKind::DOUBLE);
  EXPECT_TRUE(getControlParameterValueKindFromString("3") == ControlParameterValueKind::S64);
  EXPECT_TRUE(getControlParameterValueKindFromString("[1,2,3]") == ControlParameterValueKind::VEC3_DOUBLE);
  EXPECT_TRUE(getControlParameterValueKindFromString("[1.f,2.f,3.f]") == ControlParameterValueKind::VEC3_FLOAT);
}

TEST(ControlParams, DefineFromYaml) {
  ControlParameters params("robot-parameters");
  params.defineAndInitializeFromYamlFile(getConfigDirectoryPath() +
CHEETAH_3_DEFAULT_PARAMETERS);

  params.collection.deleteAll();
  params.defineAndInitializeFromYamlFile(getConfigDirectoryPath() + MINI_CHEETAH_DEFAULT_PARAMETERS);

  //printf("Result: \n %s\n", params.collection.printToIniString().c_str());
}