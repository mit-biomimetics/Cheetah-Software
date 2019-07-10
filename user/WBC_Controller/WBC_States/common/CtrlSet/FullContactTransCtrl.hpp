#ifndef FULL_CONTACT_TRANSITION_CTRL
#define FULL_CONTACT_TRANSITION_CTRL

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
class FullContactTransCtrl : public Controller<T> {
 public:
  FullContactTransCtrl(const FloatingBaseModel<T>*);
  virtual ~FullContactTransCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& setting_file_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  Task<T>* _body_pos_task;
  Task<T>* _body_ori_task;

  ContactSpec<T>* _fr_contact;
  ContactSpec<T>* _fl_contact;
  ContactSpec<T>* _hr_contact;
  ContactSpec<T>* _hl_contact;

  KinWBC<T>* _kin_wbc;
  WBLC<T>* _wblc;
  WBLC_ExtraData<T>* _wblc_data;

  DVec<T> _base_pos_ini;
  Vec3<T> _ini_base_pos;

  DVec<T> _Kp;
  DVec<T> _Kd;

  std::vector<T> _Kp_joint, _Kd_joint;

  DVec<T> _ini_jpos;
  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  bool _b_set_height_target;
  T _end_time;
  T _target_body_height;
  T _ini_body_height;
  Vec3<T> _ini_body_pos;

  T _max_rf_z;
  T _min_rf_z;
  int _dim_contact;
  T _ctrl_start_time;

  void _task_setup();
  void _contact_setup();
  void _compute_torque_wblc(DVec<T>& gamma);

  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
  std::string _test_file_name;
};

#endif
