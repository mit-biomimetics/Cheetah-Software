#ifndef ROBOT_CONTROLLER_H
#define ROBOT_CONTROLLER_H

#include "Controllers/LegController.h"
#include "Dynamics/FloatingBaseModel.h"
#include <Controllers/StateEstimatorContainer.h>
#include <SimUtilities/VisualizationData.h>
#include "SimUtilities/GamepadCommand.h"

class RobotController{
  friend class RobotRunner;
public:
  RobotController(){}
  virtual ~RobotController(){}

  virtual void initializeController() = 0;
/**
 * Called one time every control loop 
 */
  virtual void runController() = 0;
  virtual void updateVisualization() = 0;
  virtual ControlParameters* getUserControlParameters() = 0;

protected:
  Quadruped<float>* _quadruped;
  FloatingBaseModel<float>* _model;
  LegController<float>* _legController;
  StateEstimatorContainer<float>* _stateEstimator;
  StateEstimate<float>* _stateEstimate;
  GamepadCommand* _driverCommand;
  RobotControlParameters* _controlParameters;

  VisualizationData* _visualizationData;
  RobotType _robotType;
};

#endif
