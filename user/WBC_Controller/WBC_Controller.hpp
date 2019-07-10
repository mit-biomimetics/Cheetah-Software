#ifndef WBC_CONTROLLER
#define WBC_CONTROLLER

#include <RobotController.h>
#include "WBC_States/Cheetah_DynaCtrl_Definition.h"
#include <gui_main_control_settings_t.hpp>

template<typename T> class Test;

class WBC_Controller: public RobotController{
public:
  WBC_Controller();
  virtual ~WBC_Controller(){}
  // For WBC state (test)
  virtual void initializeController();
  virtual void runController();
  virtual void updateVisualization();
  virtual ControlParameters* getUserControlParameters() {
    return nullptr;
  }

protected:
  Test<float>* _wbc_state;
  Cheetah_Data<float>* _data;
  Cheetah_Extra_Data<float>* _extra_data;
  gui_main_control_settings_t main_control_settings;

  void _StepLocationVisualization();
  void _BodyPathVisualization();
  void _BodyPathArrowVisualization();
  void _testDebugVisualization();

  float _ini_yaw;
};



#endif
