/*!
 * @file RobotController.h
 * @brief Common framework for running robot controllers.
 * This code is a common interface between control code and hardware/simulation for mini cheetah and cheetah 3
 */

#ifndef PROJECT_ROBOTCONTROLLER_H
#define PROJECT_ROBOTCONTROLLER_H

#include <SimUtilities/IMUTypes.h>
#include <ControlParameters/ControlParameterInterface.h>
#include <ControlParameters/RobotParameters.h>
#include "Controllers/LegController.h"
#include "Dynamics/Quadruped.h"
#include "SimUtilities/GamepadCommand.h"
#include "SimUtilities/VisualizationData.h"
#include <Controllers/StateEstimatorContainer.h>
#include "cheetah_visualization_lcmt.hpp"
#include "Controllers/GaitScheduler.h"
#include "Controllers/ContactEstimator.h"
#include "Controllers/DesiredStateCommand.h"
#include "JPosInitializer.h"
#include "Utilities/PeriodicTask.h"
#include "ControlFSM.h"
#include <gui_main_control_settings_t.hpp>


/**
 * Enumerate all of the Control logic options. Each one will have 
 * an initialize and a run function that pieces together the control
 * logic from the various robot controllers and sensors.
 */
enum class ControlLogicName {
  CONTROL_FSM,
  WBC
};


#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
template <typename T> class Test;

class RobotController : public PeriodicTask {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using PeriodicTask::PeriodicTask;
  void init() override;
  void run() override;
  void cleanup() override;

  // Control options logic initializations 
  void initializeControlOptionControlFSM(); // ControlFSM
  void initializeControlOptionWBC();        // WBC

  // Control options logic run functions
  void runControlOptionControlFSM();  // ControlFSM
  void runControlOptionWBC();         // WBC

  // Initialize the state estimator with default no cheaterMode
  void initializeStateEstimator(bool cheaterMode = false);
  virtual ~RobotController();


  GamepadCommand* driverCommand;
  RobotType  robotType;
  KvhImuData* kvhImuData;
  VectorNavData* vectorNavData;
  CheaterState<double>* cheaterState;
  SpiData* spiData;
  SpiCommand* spiCommand;
  TiBoardCommand* tiBoardCommand;
  TiBoardData* tiBoardData;
  RobotControlParameters* controlParameters;
  VisualizationData* visualizationData;
  CheetahVisualization* cheetahMainVisualization;

private:
  float _ini_yaw;

  int iter = 0;

  void setupStep();
  void finalizeStep();
  void testDebugVisualization();
  void StepLocationVisualization();
  void BodyPathVisualization();
  void BodyPathArrowVisualization();

  JPosInitializer<float>* _jpos_initializer;
  Quadruped<float> _quadruped;
  LegController<float>* _legController = nullptr;
  StateEstimate<float> _stateEstimate;
  StateEstimatorContainer<float>* _stateEstimator;
  bool _cheaterModeEnabled = false;
  DesiredStateCommand<float>* _desiredStateCommand;
  ControlFSM<float>* _controlFSM;
  gui_main_control_settings_t main_control_settings;

  // Option for the control logic
  ControlLogicName CONTROL_OPTION;

  // Gait Scheduler controls the nominal contact schedule for the feet
  GaitScheduler<float>* _gaitScheduler;

  // Contact Estimator to calculate estimated forces and contacts
  ContactEstimator<double>* _contactEstimator;

  // For WBC state (test)
  Test<float>* _wbc_state;
  Cheetah_Data<float>* _data;
  Cheetah_Extra_Data<float>* _extra_data;
  FloatingBaseModel<float> _model;
  u64 _iterations = 0;
};


#endif //PROJECT_ROBOTCONTROLLER_H
