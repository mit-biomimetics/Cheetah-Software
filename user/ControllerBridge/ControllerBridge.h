#ifndef CONTROLLER_BRIDGE_H
#define CONTROLLER_BRIDGE_H

#include <RobotController.h>
#include <lcm/lcm-cpp.hpp>

#include <thread>
#include <mutex>

class ControllerBridge : public RobotController {
public:
  ControllerBridge():RobotController(), _lcm(getLcmUrl(255)) {

    _lcm.subscribe("ctrl_bridge", &ControllerBridge::handleLcm, this);
    _lcmThread = std::thread([&](){
      for(;;) _lcm.handle();
    });
  }
  virtual ~ControllerBridge(){}

  virtual void initializeController(){}
  virtual void runController();
  virtual void updateVisualization(){}
  virtual ControlParameters* getUserControlParameters() {
    return nullptr;
  }

  void handleLcm(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
    const leg_control_command_lcmt* msg) {
    (void)rbuf;
    (void)chan;
    _mutex.lock();
    command = *msg;
    _mutex.unlock();
  }

private:
  lcm::LCM _lcm;
  leg_control_command_lcmt command;
  std::thread _lcmThread;
  std::mutex _mutex;
};

#endif 