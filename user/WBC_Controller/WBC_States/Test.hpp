#ifndef TEST_H
#define TEST_H

#include "Cheetah_DynaCtrl_Definition.h"
#include "Controller.hpp"
#include "StateProvider.hpp"

template <typename T>
class Test {
 public:
  Test(FloatingBaseModel<T>* robot, const RobotType& robot_type);
  virtual ~Test();

  constexpr static T dt = 0.001;
  void GetCommand(const Cheetah_Data<T>* data, LegControllerCommand<T>* command,
                  Cheetah_Extra_Data<T>* ext_data);

  void ComputeCommand(void* _command);

  int getPhase() { return _phase; }

 protected:
  virtual void _UpdateTestOneStep() {}
  virtual void _TestInitialization() = 0;
  virtual int _NextPhase(const int& phase) = 0;
  virtual void _UpdateExtraData(Cheetah_Extra_Data<T>* ext_data) = 0;

  bool _b_first_visit;
  bool _b_save_file;
  int _phase;
  std::vector<Controller<T>*> _state_list;

  FloatingBaseModel<T>* _robot;
  RobotType _robot_type;
  FBModelState<T> _state;
  int _count;
  int _waiting_count;

  LegControllerCommand<T>* _copy_cmd;

  T _roll_limit;
  T _pitch_limit;
  Vec3<T> _body_rpy;
  bool _b_running;

  void _ParameterSetting();
  bool _Initialization(const Cheetah_Data<T>*, LegControllerCommand<T>*);

  void _SafetyCheck();
  void _SafeCommand(const Cheetah_Data<T>*, LegControllerCommand<T>*);

  StateProvider<T>* _sp;
};

template <typename T>
constexpr T Test<T>::dt;
#endif
