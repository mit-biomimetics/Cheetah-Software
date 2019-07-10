#ifndef PROJECT_LEGINVDYNSUSERPARAMETERS_H
#define PROJECT_LEGINVDYNSUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class Leg_InvDyn_UserParameters : public ControlParameters {
public:
  Leg_InvDyn_UserParameters()
      : ControlParameters("user-parameters"),
        INIT_PARAMETER(num_moving_legs)
      {}

  DECLARE_PARAMETER(double, num_moving_legs);
};

#endif //PROJECT_LEGINVDYNSUSERPARAMETERS_H
