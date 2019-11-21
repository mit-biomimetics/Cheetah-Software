/*!
 * @file RobotInterface.h
 * @brief Interface between simulator and hardware using LCM.
 */

#ifndef PROJECT_ROBOTINTERFACE_H
#define PROJECT_ROBOTINTERFACE_H

#include <ControlParameters/RobotParameters.h>
#include <Dynamics/Quadruped.h>
#include <Utilities/PeriodicTask.h>
#include <cheetah_visualization_lcmt.hpp>
#include <condition_variable>
#include <lcm-cpp.hpp>
#include <mutex>
#include <thread>
#include "Graphics3D.h"
#include "control_parameter_request_lcmt.hpp"
#include "control_parameter_respones_lcmt.hpp"
#include "gamepad_lcmt.hpp"

#define ROBOT_INTERFACE_UPDATE_PERIOD (1.f / 60.f)
#define INTERFACE_LCM_NAME "interface"
#define TIMES_TO_RESEND_CONTROL_PARAM 5

class RobotInterface : PeriodicTask {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RobotInterface(RobotType robotType, Graphics3D* gfx, PeriodicTaskManager* tm, ControlParameters& userParameters);
  RobotControlParameters& getParams() { return _controlParameters; }
  void startInterface();
  void stopInterface();
  void lcmHandler();
  void sendControlParameter(const std::string& name,
                            ControlParameterValue value,
                            ControlParameterValueKind kind, bool isUser);

  void handleControlParameter(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const control_parameter_respones_lcmt* msg);

  void handleVisualizationData(const lcm::ReceiveBuffer* rbuf,
                               const std::string& chan,
                               const cheetah_visualization_lcmt* msg);

  void init() {}
  void run();
  void cleanup() {}
  virtual ~RobotInterface() {
    delete _simulator;
    stop();
  }

 private:
  PeriodicTaskManager _taskManager;
  gamepad_lcmt _gamepad_lcmt;
  control_parameter_request_lcmt _parameter_request_lcmt;
  bool _pendingControlParameterSend = false;
  lcm::LCM _lcm;
  uint64_t _robotID;
  std::thread _lcmThread;
  VisualizationData _visualizationData;
  RobotControlParameters _controlParameters;
  ControlParameters& _userParameters;
  Graphics3D* _gfx;
  RobotType _robotType;
  bool _running = false;

  std::mutex _lcmMutex;
  std::condition_variable _lcmCV;
  bool _waitingForLcmResponse = false;
  bool _lcmResponseBad = true;

  // forward kinematics
  Quadruped<double> _quadruped;
  FloatingBaseModel<double> _model;
  DynamicsSimulator<double>* _simulator = nullptr;
  FBModelState<double> _fwdKinState;
};

#endif  // PROJECT_ROBOTINTERFACE_H
