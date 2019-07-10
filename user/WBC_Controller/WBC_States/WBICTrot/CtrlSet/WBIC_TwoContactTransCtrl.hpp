#ifndef WBIC_TROT_TWO_CONTACT_TRANSITION_CHEETAH
#define WBIC_TROT_TWO_CONTACT_TRANSITION_CHEETAH

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
class WBIC_TwoContactTransCtrl : public Controller<T> {
 public:
  WBIC_TwoContactTransCtrl(WBICTrotTest<T>* test,
                           const FloatingBaseModel<T>* robot, size_t cp1,
                           size_t cp2, int transit_dir);
  virtual ~WBIC_TwoContactTransCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  virtual void CtrlInitialization(const std::string& setting_file_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  void _SetContact(const size_t& cp_idx, const T& upper_lim);

  WBICTrotTest<T>* _trot_test;
  size_t _cp1, _cp2;
  int _transit_dir;  // 1: swing start, -1: swing end

  Task<T>* _body_pos_task;
  Task<T>* _body_ori_task;

  ContactSpec<T>* fr_contact_;
  ContactSpec<T>* fl_contact_;
  ContactSpec<T>* hr_contact_;
  ContactSpec<T>* hl_contact_;

  KinWBC<T>* kin_wbc_;
  WBIC<T>* wbic_;
  WBIC_ExtraData<T>* wbic_data_;

  DVec<T> base_pos_ini_;
  Vec3<T> ini_base_pos_;

  std::vector<T> _Kp_joint, _Kd_joint;

  DVec<T> ini_jpos_;

  bool b_set_height_target_;
  T end_time_;
  T _body_height_cmd;
  T ini_body_height_;
  Vec3<T> ini_body_pos_;

  T max_rf_z_;
  T min_rf_z_;
  int _dim_contact;
  T ctrl_start_time_;

  void _task_setup();
  void _contact_setup();
  void _compute_torque_wbic(DVec<T>& gamma);

  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
  std::string _test_file_name;
};
#endif
