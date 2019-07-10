#ifndef BACKFLIP_CTRL
#define BACKFLIP_CTRL

#include <ParamHandler/ParamHandler.hpp>
#include <WBC_States/BackFlip/DataReader.hpp>
#include <WBC_States/Controller.hpp>

template <typename T>
class ContactSpec;
template <typename T>
class WBLC;
template <typename T>
class WBLC_ExtraData;
template <typename T>
class KinWBC;
template <typename T>
class StateProvider;
template <typename T>
class WBLCTrotTest;

class DataReader;

template <typename T>
class BackFlipCtrl : public Controller<T> {
 public:
  BackFlipCtrl(const FloatingBaseModel<T>*, DataReader*);
  virtual ~BackFlipCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  DataReader* _data_reader;

  DVec<T> _Kp, _Kd;
  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _jtorque;

  std::vector<T> _Kp_joint, _Kd_joint;

  bool _b_set_height_target;
  T _end_time;
  int _dim_contact;

  void _update_joint_command();

  T _ctrl_start_time;
  T _q_knee_max;
  T _qdot_knee_max;

  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
};

#endif
