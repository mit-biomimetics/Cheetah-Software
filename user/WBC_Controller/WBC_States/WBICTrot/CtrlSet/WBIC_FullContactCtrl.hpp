#ifndef WBIC_TROT_FULL_CONTACT_BODY_CTRL
#define WBIC_TROT_FULL_CONTACT_BODY_CTRL

#include <ParamHandler/ParamHandler.hpp>
#include <WBC_States/Controller.hpp>

template <typename T>
class ContactSpec;
template <typename T>
class WBIC;
template <typename T>
class WBIC_ExtraData;
template <typename T>
class KinWBC;
template <typename T>
class StateProvider;
template <typename T>
class WBICTrotTest;

template <typename T>
class WBIC_FullContactCtrl : public Controller<T> {
 public:
  WBIC_FullContactCtrl(WBICTrotTest<T>*, const FloatingBaseModel<T>*);
  virtual ~WBIC_FullContactCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  WBICTrotTest<T>* _trot_test;
  DVec<T> _Kp, _Kd;
  std::vector<T> _Kp_joint;
  std::vector<T> _Kd_joint;

  bool _b_set_height_target;
  T _end_time;
  int _dim_contact;

  Task<T>* _body_pos_task;
  Task<T>* _body_ori_task;

  KinWBC<T>* kin_wbc_;
  ContactSpec<T>* fr_contact_;
  ContactSpec<T>* fl_contact_;
  ContactSpec<T>* hr_contact_;
  ContactSpec<T>* hl_contact_;
  WBIC<T>* wbic_;
  WBIC_ExtraData<T>* wbic_data_;

  T target_body_height_;
  Vec3<T> _ini_rpy;
  T ini_body_height_;

  void _task_setup();
  void _contact_setup();
  void _compute_torque_wbic(DVec<T>& gamma);

  T _ctrl_start_time;
  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
};

#endif
