#ifndef BOUNDING_CTRL
#define BOUNDING_CTRL

#include <WBC_Ctrl/WBC_Ctrl.hpp>
#include "ImpulseCurve.hpp"

class MIT_UserParameters;

template <typename T>
class BoundingCtrl : public WBC_Ctrl<T> {
 public:
  BoundingCtrl(FloatingBaseModel<T> );
  virtual ~BoundingCtrl();

 protected:
  virtual void _ContactTaskUpdate(void * input, ControlFSMData<T> & data);
  void _body_task_setup();
  void _leg_task_setup();
  void _contact_update();
  void _FirstVisit();

  void _StatusCheck();
  void _setupTaskAndContactList();
  void _ContactUpdate();

  void _lcm_data_sending();
  void _save_file();

  T _curr_time;
  T dt;

  bool _b_front_swing;
  bool _b_hind_swing;

  bool _b_front_contact_est;
  bool _b_hind_contact_est;

  T _step_width;

  T _contact_vel_threshold;
  T _K_time = 0.0;
  T _swing_height;

  T _front_foot_offset;
  T _hind_foot_offset;
 
  T _aerial_duration;
  T _front_swing_time;

  T _front_previous_stance;
  T _front_previous_swing;
  T _front_current_stance;

  T _hind_previous_stance;
  T _hind_previous_swing;
  T _hind_current_stance;

  T _nominal_gait_period;
  T _gait_period;

  T _swing_time = 0.18;
  T _stance_time;
  T _nominal_stance_time = 0.15;
  T _step_length_lim = 0.3;

  int _dim_contact;

  T _front_start_time;
  T _hind_start_time;

  T _front_time;
  T _hind_time;

  Task<T>* _local_head_pos_task;
  Task<T>* _local_tail_pos_task;

  Task<T>* _local_roll_task;
  Task<T>* _body_ryrz_task;

  Task<T>* _fr_foot_local_task;
  Task<T>* _fl_foot_local_task;
  Task<T>* _hr_foot_local_task;
  Task<T>* _hl_foot_local_task;

  Vec3<T> _fr_foot_pos;
  DVec<T> _fr_foot_vel;
  DVec<T> _fr_foot_acc;

  Vec3<T> _fl_foot_pos;
  DVec<T> _fl_foot_vel;
  DVec<T> _fl_foot_acc;

  Vec3<T> _hr_foot_pos;
  DVec<T> _hr_foot_vel;
  DVec<T> _hr_foot_acc;

  Vec3<T> _hl_foot_pos;
  DVec<T> _hl_foot_vel;
  DVec<T> _hl_foot_acc;

  ContactSpec<T>* _fr_contact;
  ContactSpec<T>* _fl_contact;
  ContactSpec<T>* _hr_contact;
  ContactSpec<T>* _hl_contact;

  T _target_leg_height = 0.28;
  T _total_mass;

  ImpulseCurve<T> _front_z_impulse;
  ImpulseCurve<T> _hind_z_impulse;

  Vec3<T> _ini_front_body;
  Vec3<T> _ini_hind_body;

  Vec3<T> _ini_fr;
  Vec3<T> _ini_fl;
  Vec3<T> _ini_hr;
  Vec3<T> _ini_hl;

  Vec3<T> _fin_fr;
  Vec3<T> _fin_fl;
  Vec3<T> _fin_hr;
  Vec3<T> _fin_hl;

  Vec3<T> _head_pos_ini;
  Vec3<T> _tail_pos_ini;

  Vec3<T> _vel_des;
  DVec<T> _Fr_des;
};

#endif
