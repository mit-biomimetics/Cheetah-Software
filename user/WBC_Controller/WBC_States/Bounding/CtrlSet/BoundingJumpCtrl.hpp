#ifndef BOUNDING_JUMP_CTRL
#define BOUNDING_JUMP_CTRL

#include <WBC_States/Controller.hpp>
#include <casadi/casadi.hpp>
#include <ParamHandler/ParamHandler.hpp>

template <typename T> class ContactSpec;
template <typename T> class KinWBC;
template <typename T> class WBIC;
template <typename T> class WBIC_ExtraData;
template <typename T> class StateProvider;
template <typename T> class BoundingTest;

using namespace casadi;

template <typename T>
class BoundingJumpCtrl: public Controller<T> {
 public:
  BoundingJumpCtrl(BoundingTest<T>*, const FloatingBaseModel<T>*);
  virtual ~BoundingJumpCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  /*!
   * Parameter setting from yaml file
   * @param category_name : name of this controller
   */
  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  void _update_motion_cmd();
  int _iter;
  int _motion_count;
  MX simple_dyn(const MX& x, const MX& u, const MX& P, int i);

  void _optimize_jump(const DVec<T> & X0, T xdot_target, T zdot_target );
  void _contact_update();
  void _compute_torque_wbic(DVec<T> & gamma);

  std::vector<T> _Kp_joint, _Kd_joint;

  KinWBC<T>* _kin_wbc;
  WBIC<T>* _wbic;
  WBIC_ExtraData<T>* _wbic_data;

  T _ctrl_start_time;

  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  BoundingTest<T>* _bounding_test;
  T _target_leg_height;

  std::vector<T> _xpos_list;
  std::vector<T> _zpos_list;
  std::vector<T> _theta_list;

  std::vector<T> _xvel_list;
  std::vector<T> _zvel_list;
  std::vector<T> _theta_dot_list;

  std::vector<T> _fr_x_list; 
  std::vector<T> _fr_z_list; 

  std::vector<T> _P_foot;

  Task<T>* _local_roll_task;
  Task<T>* _body_ryrz_task;

  Task<T>* _fr_foot_local_task;
  Task<T>* _fl_foot_local_task;
  Task<T>* _hr_foot_local_task;
  Task<T>* _hl_foot_local_task;

  Task<T>* _local_head_pos_task;
  Task<T>* _local_tail_pos_task;

  ContactSpec<T>* _fr_contact;
  ContactSpec<T>* _fl_contact; 
  ContactSpec<T>* _hr_contact; 
  ContactSpec<T>* _hl_contact; 



  T _z_vel_target;
  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
};
#endif
