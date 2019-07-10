#ifndef MIT_CONTROLLER
#define MIT_CONTROLLER

#include <RobotController.h>
#include "Controllers/GaitScheduler.h"
#include "Controllers/ContactEstimator.h"
#include "Controllers/DesiredStateCommand.h"
#include "FSM_States/ControlFSM.h"
#include "MIT_UserParameters.h"
#include <gui_main_control_settings_t.hpp>

class MIT_Controller: public RobotController{
public:
  MIT_Controller();
  virtual ~MIT_Controller(){}

  virtual void initializeController();
  virtual void runController();
  virtual void updateVisualization(){}
  virtual ControlParameters* getUserControlParameters() {
    return &userParameters;
  }


protected:
  ControlFSM<float>* _controlFSM;
  // Gait Scheduler controls the nominal contact schedule for the feet
  DesiredStateCommand<float>* _desiredStateCommand;
  GaitScheduler<float>* _gaitScheduler;
  MIT_UserParameters userParameters;

  gui_main_control_settings_t main_control_settings;
};


#endif
