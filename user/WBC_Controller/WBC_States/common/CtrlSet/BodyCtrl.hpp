#ifndef BODY_CTRL
#define BODY_CTRL

#include <ParamHandler/ParamHandler.hpp>
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
class BodyCtrlTest;

template <typename T>
class BodyCtrl : public Controller<T> {
 public:
  BodyCtrl(const FloatingBaseModel<T>*);
  virtual ~BodyCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

  Vec3<T> _target_ori_command;

 protected:
  BodyCtrlTest<T>* _body_test;
  std::vector<T> _Kp_joint, _Kd_joint;
  DVec<T> _Kp, _Kd;
  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  bool _b_set_height_target;
  T _end_time;
  int _dim_contact;

  Task<T>* _body_pos_task;
  Task<T>* _body_ori_task;

  KinWBC<T>* _kin_wbc;
  ContactSpec<T>* _fr_contact;
  ContactSpec<T>* _fl_contact;
  ContactSpec<T>* _hr_contact;
  ContactSpec<T>* _hl_contact;
  WBLC<T>* _wblc;
  WBLC_ExtraData<T>* _wblc_data;

  T _target_body_height;

  void _task_setup();
  void _contact_setup();
  void _compute_torque_wblc(DVec<T>& gamma);

  T _ctrl_start_time;
  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
};

#endif
