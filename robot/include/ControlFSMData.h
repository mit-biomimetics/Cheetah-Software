#ifndef CONTROLFSMDATA_H
#define CONTROLFSMDATA_H

#include "Dynamics/Quadruped.h"
#include "Controllers/StateEstimatorContainer.h"
#include "Controllers/LegController.h"
#include "Controllers/GaitScheduler.h"
#include "Controllers/DesiredStateCommand.h"
#include <ControlParameters/RobotParameters.h>


/**
 *
 */
template <typename T>
struct ControlFSMData {
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Quadruped<T>* _quadruped;
  StateEstimatorContainer<T>* _stateEstimator;
  LegController<T>* _legController;
  GaitScheduler<T>* _gaitScheduler;
  DesiredStateCommand<T>* _desiredStateCommand;
  RobotControlParameters* controlParameters;
};


template struct ControlFSMData<double>;
template struct ControlFSMData<float>;


#endif // CONTROLFSM_H