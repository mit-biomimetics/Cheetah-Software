#ifndef POSTURE_KEEPING_CTRL
#define POSTURE_KEEPING_CTRL

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
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
class PostureKeepingCtrl : public Controller<T> {
 public:
  PostureKeepingCtrl(const FloatingBaseModel<T>*);
  virtual ~PostureKeepingCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

  void SetJPosTarget(const DVec<T>& jpos_target) { _jpos_target = jpos_target; }

 protected:
  DVec<T> _Kp, _Kd;
  std::vector<T> _Kp_joint, _Kd_joint;

  bool _b_target_set;

  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  DVec<T> _jpos_ini;
  DVec<T> _jpos_target;

  T _end_time;
  int _dim_contact;
  T _ctrl_start_time;

  ContactSpec<T>* _contact;
  WBLC<T>* _wblc;
  WBLC_ExtraData<T>* _wblc_data;

  void _task_setup();
  void _contact_setup();
  void _compute_torque_wblc(DVec<T>& gamma);

  ParamHandler* _param_handler = NULL;
  StateProvider<T>* _sp;
};

#endif
