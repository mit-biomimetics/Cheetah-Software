#ifndef JPOS_CONTROLLER
#define JPOS_CONTROLLER

#include <RobotController.h>
#include "Leg_InvDyn_UserParameters.h"

class Leg_InvDyn_Controller:public RobotController{
  public:
    Leg_InvDyn_Controller():RobotController(){
    }
    virtual ~Leg_InvDyn_Controller(){}

    virtual void initializeController(){}
    virtual void runController();
    virtual void updateVisualization(){}
    virtual ControlParameters* getUserControlParameters() {
      return &userParameters;
    }
  protected:
    Leg_InvDyn_UserParameters userParameters;
};

#endif
