#ifndef PROJECT_JPOSUSERPARAMETERS_H
#define PROJECT_JPOSUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class JPosUserParameters : public ControlParameters {
public:
  JPosUserParameters()
      : ControlParameters("user-parameters"),
        INIT_PARAMETER(testValue),
        INIT_PARAMETER(testValue2)
      {}

  DECLARE_PARAMETER(double, testValue);
  DECLARE_PARAMETER(double, testValue2)
};

#endif //PROJECT_JPOSUSERPARAMETERS_H
