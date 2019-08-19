#ifndef PROJECT_MITUSERPARAMETERS_H
#define PROJECT_MITUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class MIT_UserParameters : public ControlParameters {
public:
  MIT_UserParameters()
      : ControlParameters("user-parameters"),
        INIT_PARAMETER(cmpc_gait),
        INIT_PARAMETER(test_param)
      {}

  DECLARE_PARAMETER(double, cmpc_gait);
  DECLARE_PARAMETER(double, test_param);
};

#endif //PROJECT_MITUSERPARAMETERS_H
