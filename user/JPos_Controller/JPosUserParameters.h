#ifndef PROJECT_JPOSUSERPARAMETERS_H
#define PROJECT_JPOSUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class JPosUserParameters : public ControlParameters {
public:
  JPosUserParameters()
      : ControlParameters("jpos-parameters"),
        INIT_PARAMETER(cmpc_gait),
        INIT_PARAMETER(test_param){
        }

  DECLARE_PARAMETER(double, cmpc_gait);
  DECLARE_PARAMETER(double, test_param);
};

#endif //PROJECT_JPOSUSERPARAMETERS_H
