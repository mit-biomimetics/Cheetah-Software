#include "BoundingCtrl.hpp"

#include <WBC_Ctrl/TaskSet/BodyRyRzTask.hpp>
#include <WBC_Ctrl/TaskSet/LocalHeadPosTask.hpp>
#include <WBC_Ctrl/TaskSet/LocalPosTask.hpp>
#include <WBC_Ctrl/TaskSet/LocalRollTask.hpp>
#include <WBC_Ctrl/TaskSet/LocalTailPosTask.hpp>

#include <WBC_Ctrl/ContactSet/SingleContact.hpp>

template <typename T>
BoundingCtrl<T>::BoundingCtrl(FloatingBaseModel<T> robot):
  WBC_Ctrl<T>(robot),
  _curr_time(0.),
      _step_width(0.07),
      _K_time(0.5),
      _swing_height(0.08),
      _front_foot_offset(0.02),
      _hind_foot_offset(-0.05),
      _dim_contact(0),
      _fr_foot_vel(3),
      _fr_foot_acc(3),
      _fl_foot_vel(3),
      _fl_foot_acc(3),
      _hr_foot_vel(3),
      _hr_foot_acc(3),
      _hl_foot_vel(3),
      _hl_foot_acc(3),
      _total_mass(9.){

  _vel_des.setZero();
  
  // Start from front swing & hind stance
  _b_front_swing = true;
  _b_hind_swing = false;

  _fr_foot_vel.setZero();
  _fr_foot_acc.setZero();

  _fl_foot_vel.setZero();
  _fl_foot_acc.setZero();

  _hr_foot_vel.setZero();
  _hr_foot_acc.setZero();

  _hl_foot_vel.setZero();
  _hl_foot_acc.setZero();

  _local_roll_task = new LocalRollTask<T>(&(WBCtrl::_model));
  _body_ryrz_task = new BodyRyRzTask<T>(&(WBCtrl::_model));

  _fr_foot_local_task = new LocalPosTask<T>(&(WBCtrl::_model), linkID::FR, linkID::FR_abd);
  _fl_foot_local_task = new LocalPosTask<T>(&(WBCtrl::_model), linkID::FL, linkID::FL_abd);
  _hr_foot_local_task = new LocalPosTask<T>(&(WBCtrl::_model), linkID::HR, linkID::HR_abd);
  _hl_foot_local_task = new LocalPosTask<T>(&(WBCtrl::_model), linkID::HL, linkID::HL_abd);

  _local_head_pos_task = new LocalHeadPosTask<T>(&(WBCtrl::_model));
  _local_tail_pos_task = new LocalTailPosTask<T>(&(WBCtrl::_model));

  //((LocalHeadPosTask<T>*)_local_head_pos_task)->_Kp[0] = 0.;
  //((LocalHeadPosTask<T>*)_local_head_pos_task)->_Kp[1] = 0.;

  //((LocalTailPosTask<T>*)_local_tail_pos_task)->_Kp[0] = 0.;
  //((LocalTailPosTask<T>*)_local_tail_pos_task)->_Kp[1] = 0.;

  _fr_contact = new SingleContact<T>(&(WBCtrl::_model), linkID::FR);
  _fl_contact = new SingleContact<T>(&(WBCtrl::_model), linkID::FL);
  _hr_contact = new SingleContact<T>(&(WBCtrl::_model), linkID::HR);
  _hl_contact = new SingleContact<T>(&(WBCtrl::_model), linkID::HL);
}

template <typename T>
void BoundingCtrl<T>::_ContactUpdate() {
  if (_b_front_swing) {  _b_front_contact_est = false;
  } else {  _b_front_contact_est = true; }

  if (_b_hind_swing) { _b_hind_contact_est = false;
  } else {  _b_hind_contact_est = true;  }

  if (_b_front_swing && (_front_time > 0.6 * _swing_time)) {
    if (_front_time > _swing_time) {
      _b_front_contact_est = true;
    }
  }
  if (_b_hind_swing && (_hind_time > 0.6 * _swing_time)) {
    if (_hind_time > _swing_time) {
      _b_hind_contact_est = true;
    }
  }
}

template <typename T>
void BoundingCtrl<T>::_StatusCheck() {
  _stance_time = _nominal_stance_time;
  _gait_period = _nominal_gait_period;

  // High speed bounding stance time scaling
  if (fabs(_vel_des[0]) * _nominal_stance_time >
      _step_length_lim) {
    _stance_time = fabs(_step_length_lim / _vel_des[0]);
    _gait_period = _stance_time + _swing_time;
  }

  T adjust_time(0.);
  if (_b_front_swing && (_front_time > 0.5 * _swing_time)) {
    // Check Contact
    if (_b_front_contact_est) {
      _b_front_swing = false;
      _front_start_time = _curr_time;

      _front_previous_swing = _front_time;
      adjust_time = _K_time * (_hind_previous_stance + _aerial_duration -
                               0.5 * _gait_period);
      adjust_time = coerce<T>(adjust_time, 0, 0.04);
      _front_current_stance = _stance_time - adjust_time;

      T apex = _total_mass * 9.81 * (_swing_time + _front_current_stance) /
               (2. * 2.0 * 0.7 * _front_current_stance);


      _front_z_impulse.setCurve(apex, _front_current_stance);
      _front_previous_stance = _front_current_stance;

      _ini_front_body = 0.5 * (WBCtrl::_model)._pGC[linkID::FR_abd] +
                        0.5 * (WBCtrl::_model)._pGC[linkID::FL_abd] -
                        0.5 * (WBCtrl::_model)._pGC[linkID::FR] -
                        0.5 * (WBCtrl::_model)._pGC[linkID::FL];

      _front_swing_time = _swing_time - 2. * dt;

      // Front time reset
      _front_time = 0.;
    }
  }

  if (_b_hind_swing && (_hind_time > 0.5 * _swing_time)) {
    // Check Contact
    if (_b_hind_contact_est) {
      _b_hind_swing = false;
      _hind_start_time = _curr_time;

      _hind_previous_swing = _hind_time;
      adjust_time = _K_time * (_front_previous_stance + _aerial_duration -
                               0.5 * _gait_period);
      adjust_time = coerce<T>(adjust_time, 0., 0.04);
      _hind_current_stance = _stance_time - adjust_time;

      T apex = _total_mass * 9.81 * (_swing_time + _hind_current_stance) /
               (2. * 2.0 * 0.7 * _hind_current_stance);


      _hind_z_impulse.setCurve(apex, _hind_current_stance);
      _hind_previous_stance = _hind_current_stance;

      _ini_hind_body = 0.5 * (WBCtrl::_model)._pGC[linkID::HR_abd] +
                       0.5 * (WBCtrl::_model)._pGC[linkID::HL_abd] -
                       0.5 * (WBCtrl::_model)._pGC[linkID::HR] -
                       0.5 * (WBCtrl::_model)._pGC[linkID::HL];

      // Hind time reset
      _hind_time = 0.;
    }
  }

  T scale(1.0);
  // If stance time is over switch to swing
  if ((!_b_front_swing) && (_front_time > (_front_current_stance))) {
    _b_front_swing = true;

    _ini_fr = (WBCtrl::_model)._pGC[linkID::FR] -
              (WBCtrl::_model)._pGC[linkID::FR_abd];
    _ini_fl = (WBCtrl::_model)._pGC[linkID::FL] -
              (WBCtrl::_model)._pGC[linkID::FL_abd];

    _fin_fr = _vel_des * _stance_time / 2. * scale;
    _fin_fl = _vel_des * _stance_time / 2. * scale;

    _fin_fr[0] += _front_foot_offset;
    _fin_fl[0] += _front_foot_offset;

    _front_start_time = _curr_time;

    _head_pos_ini = -0.5*_ini_fr - 0.5 * _ini_fl;
    _front_time = 0.;
  }
  // If stance time is over switch to swing
  if ((!_b_hind_swing) && (_hind_time > (_hind_current_stance))) {
    _b_hind_swing = true;

    _ini_hr = (WBCtrl::_model)._pGC[linkID::HR] -
              (WBCtrl::_model)._pGC[linkID::HR_abd];
    _ini_hl = (WBCtrl::_model)._pGC[linkID::HL] -
              (WBCtrl::_model)._pGC[linkID::HL_abd];

    _fin_hr = _vel_des * _stance_time / 2. * scale;
    _fin_hl = _vel_des * _stance_time / 2. * scale;

    _fin_hr[0] += _hind_foot_offset;
    _fin_hl[0] += _hind_foot_offset;

    _hind_start_time = _curr_time;
    _tail_pos_ini = -0.5*_ini_hr - 0.5 * _ini_hl;
    _hind_time = 0.;
  }
}

template <typename T>
void BoundingCtrl<T>::_setupTaskAndContactList() {
  if (_b_front_swing) {
    WBCtrl::_task_list.push_back(_fl_foot_local_task);
    WBCtrl::_task_list.push_back(_fr_foot_local_task);
  } else {
    WBCtrl::_task_list.push_back(_local_head_pos_task);

    // Contact
    WBCtrl::_contact_list.push_back(_fr_contact);
    WBCtrl::_contact_list.push_back(_fl_contact);
  }

  if (_b_hind_swing) {
    WBCtrl::_task_list.push_back(_hr_foot_local_task);
    WBCtrl::_task_list.push_back(_hl_foot_local_task);
  } else {
    WBCtrl::_task_list.push_back(_local_tail_pos_task);

    // Contact
    WBCtrl::_contact_list.push_back(_hr_contact);
    WBCtrl::_contact_list.push_back(_hl_contact);
  }
}

template <typename T>
void BoundingCtrl<T>::_ContactTaskUpdate(void * input, ControlFSMData<T> & data) {
  (void)input;

  dt = data.controlParameters->controller_dt;
  _curr_time += dt;
  _vel_des[0] = data._desiredStateCommand->data.stateDes(6) ;
  _vel_des[1] = data._desiredStateCommand->data.stateDes(7) ;
  _vel_des[2] = 0.;

  static bool first_visit(true);
  if(first_visit){
    _curr_time = 0.;
    _FirstVisit();
    first_visit = false;
  }
  
  // Clear all
  WBCtrl::_contact_list.clear();
  WBCtrl::_task_list.clear();

  _front_time = _curr_time - _front_start_time;
  _hind_time = _curr_time - _hind_start_time;

  // Update Contact
  _ContactUpdate();
  // Update Current Stance and Swing Status
  _StatusCheck();

  if ((!_b_front_swing) || (!_b_hind_swing)) {
    _aerial_duration = 0.;
    WBCtrl::_task_list.push_back(_local_roll_task);
    WBCtrl::_task_list.push_back(_body_ryrz_task);
  } else {
    _aerial_duration += dt;
  }

  _setupTaskAndContactList();

  _contact_update();
  _body_task_setup();
  _leg_task_setup();
}

template <typename T>
void BoundingCtrl<T>::_body_task_setup() {
  // Body Ry Rz (pitch and yaw)
  Vec3<T> rpy_des;
  rpy_des.setZero();
  Quat<T> quat_des = ori::rpyToQuat(rpy_des);
  DVec<T> vel_des(2);
  vel_des.setZero();
  DVec<T> acc_des(2);
  acc_des.setZero();
  _body_ryrz_task->UpdateTask(&(quat_des), vel_des, acc_des);

  // Local Roll
  T roll_cmd(0.);
  vel_des.resize(1);
  vel_des.setZero();
  acc_des.resize(1);
  acc_des.setZero();
  _local_roll_task->UpdateTask(&(roll_cmd), vel_des, acc_des);

  // Local Head
  Vec3<T> pos_des;
  pos_des.setZero();
  pos_des[2] = _target_leg_height;
  vel_des.resize(3);
  vel_des.setZero();
  acc_des.resize(3);
  acc_des.setZero();

  if (!_b_front_swing) { 
    pos_des[0] = _head_pos_ini[0] + _vel_des[0] * _front_time;
    pos_des[1] = _head_pos_ini[1] + _vel_des[1] * _front_time;
    _local_head_pos_task->UpdateTask(&pos_des, vel_des, acc_des);
  }
  if (!_b_hind_swing) {
    pos_des[0] = _tail_pos_ini[0] + _vel_des[0] * _hind_time;
    pos_des[1] = _tail_pos_ini[1] + _vel_des[1] * _hind_time;
    _local_tail_pos_task->UpdateTask(&pos_des, vel_des, acc_des);
  }
}

template <typename T>
void BoundingCtrl<T>::_leg_task_setup() {
  // for Z (height)
  T amp(_swing_height / 2.);
  T omega(2. * M_PI / _swing_time);
  T t;

  _fr_foot_pos.setZero();
  _fl_foot_pos.setZero();
  _hr_foot_pos.setZero();
  _hl_foot_pos.setZero();

  DVec<T> vel_des(3); vel_des.setZero();
  DVec<T> acc_des(3); acc_des.setZero();

  DVec<T> Fr_des(3);
  Fr_des.setZero();

  if (_b_front_swing) {
    t = _front_time;
    if (_front_time > _swing_time - 2. * dt) {
      t = 0.;
    }

    // FR X
    _fr_foot_pos[0] = smooth_change(_ini_fr[0], _fin_fr[0], _swing_time, t);
    _fr_foot_vel[0] = smooth_change_vel(_ini_fr[0], _fin_fr[0], _swing_time, t);
    _fr_foot_acc[0] = smooth_change_acc(_ini_fr[0], _fin_fr[0], _swing_time, t);
    // FR Y
    _fr_foot_pos[1] = smooth_change(_ini_fr[1], -_step_width, _swing_time, t);
    _fr_foot_vel[1] =
        smooth_change_vel(_ini_fr[1], -_step_width, _swing_time, t);
    _fr_foot_acc[1] =
        smooth_change_acc(_ini_fr[1], -_step_width, _swing_time, t);
    // FR Z
    _fr_foot_pos[2] = -_target_leg_height + amp * (1. - cos(omega * t));
    _fr_foot_vel[2] = amp * omega * sin(omega * t);
    _fr_foot_acc[2] = amp * omega * omega * sin(omega * t);

    // FL X
    _fl_foot_pos[0] = smooth_change(_ini_fl[0], _fin_fl[0], _swing_time, t);
    _fl_foot_vel[0] = smooth_change_vel(_ini_fl[0], _fin_fl[0], _swing_time, t);
    _fl_foot_acc[0] = smooth_change_acc(_ini_fl[0], _fin_fl[0], _swing_time, t);
    // FL Y
    _fl_foot_pos[1] = smooth_change(_ini_fl[1], _step_width, _swing_time, t);
    _fl_foot_vel[1] =
        smooth_change_vel(_ini_fl[1], _step_width, _swing_time, t);
    _fl_foot_acc[1] =
        smooth_change_acc(_ini_fl[1], _step_width, _swing_time, t);
    // FL Z
    _fl_foot_pos[2] = -_target_leg_height + amp * (1. - cos(omega * t));
    _fl_foot_vel[2] = amp * omega * sin(omega * t);
    _fl_foot_acc[2] = amp * omega * omega * sin(omega * t);

    _fr_foot_local_task->UpdateTask(&(_fr_foot_pos), _fr_foot_vel, _fr_foot_acc);
    _fl_foot_local_task->UpdateTask(&(_fl_foot_pos), _fl_foot_vel, _fl_foot_acc);
  } else {  // Front Leg Stance

    WBCtrl::_wbic_data->_W_rf = DVec<T>::Constant(6, 1.0);

    WBCtrl::_wbic_data->_W_rf[2] = 50.0;
    WBCtrl::_wbic_data->_W_rf[5] = 50.0;

    Fr_des[2] = _front_z_impulse.getValue(_front_time);
    _fr_contact->setRFDesired(Fr_des);
    _fl_contact->setRFDesired(Fr_des);
  }

  if (_b_hind_swing) {
    t = _hind_time;
    if (_hind_time > _swing_time - 2. * dt) {
      t = 0.;
    }

    // HR X
    _hr_foot_pos[0] = smooth_change(_ini_hr[0], _fin_hr[0], _swing_time, t);
    _hr_foot_vel[0] = smooth_change_vel(_ini_hr[0], _fin_hr[0], _swing_time, t);
    _hr_foot_acc[0] = smooth_change_acc(_ini_hr[0], _fin_hr[0], _swing_time, t);
    // HR Y
    _hr_foot_pos[1] = smooth_change(_ini_hr[1], -_step_width, _swing_time, t);
    _hr_foot_vel[1] =
        smooth_change_vel(_ini_hr[1], -_step_width, _swing_time, t);
    _hr_foot_acc[1] =
        smooth_change_acc(_ini_hr[1], -_step_width, _swing_time, t);
    // HR Z
    _hr_foot_pos[2] = -_target_leg_height + amp * (1. - cos(omega * t));
    _hr_foot_vel[2] = amp * omega * sin(omega * t);
    _hr_foot_acc[2] = amp * omega * omega * sin(omega * t);

    // HL X
    _hl_foot_pos[0] = smooth_change(_ini_hl[0], _fin_hl[0], _swing_time, t);
    _hl_foot_vel[0] = smooth_change_vel(_ini_hl[0], _fin_hl[0], _swing_time, t);
    _hl_foot_acc[0] = smooth_change_acc(_ini_hl[0], _fin_hl[0], _swing_time, t);
    // HL Y
    _hl_foot_pos[1] = smooth_change(_ini_hl[1], _step_width, _swing_time, t);
    _hl_foot_vel[1] =
        smooth_change_vel(_ini_hl[1], _step_width, _swing_time, t);
    _hl_foot_acc[1] =
        smooth_change_acc(_ini_hl[1], _step_width, _swing_time, t);
    // HL Z
    _hl_foot_pos[2] = -_target_leg_height + amp * (1. - cos(omega * t));
    _hl_foot_vel[2] = amp * omega * sin(omega * t);
    _hl_foot_acc[2] = amp * omega * omega * sin(omega * t);

    _hr_foot_local_task->UpdateTask(&(_hr_foot_pos), _hr_foot_vel,
                                    _hr_foot_acc);
    _hl_foot_local_task->UpdateTask(&(_hl_foot_pos), _hl_foot_vel,
                                    _hl_foot_acc);

  } else {
    if (_b_front_swing) {  // Hind stance only
      WBCtrl::_wbic_data->_W_rf = DVec<T>::Constant(6, 1.0);
      WBCtrl::_wbic_data->_W_rf[2] = 50.0;
      WBCtrl::_wbic_data->_W_rf[5] = 50.0;

      Fr_des[2] = _hind_z_impulse.getValue(_hind_time);
      _hr_contact->setRFDesired(Fr_des);
      _hl_contact->setRFDesired(Fr_des);
    } else {  // Front and Hind stance
      WBCtrl::_wbic_data->_W_rf = DVec<T>::Constant(12, 1.0);
      WBCtrl::_wbic_data->_W_rf[2] = 50.0;
      WBCtrl::_wbic_data->_W_rf[5] = 50.0;
      WBCtrl::_wbic_data->_W_rf[8] = 50.0;
      WBCtrl::_wbic_data->_W_rf[11] = 50.0;

      Fr_des[2] = _front_z_impulse.getValue(_front_time);
      _fr_contact->setRFDesired(Fr_des);
      _fl_contact->setRFDesired(Fr_des);

      Fr_des[2] = _hind_z_impulse.getValue(_hind_time);
      _hr_contact->setRFDesired(Fr_des);
      _hl_contact->setRFDesired(Fr_des);
    }
  }
}

template <typename T>
void BoundingCtrl<T>::_FirstVisit() {
  T apex = _total_mass * 9.81 * 
    (_swing_time + _nominal_stance_time) /
           (2. * 2.0 * 0.7 * _nominal_stance_time);

  _front_z_impulse.setCurve(apex, _nominal_stance_time);
  _hind_z_impulse.setCurve(apex, _nominal_stance_time);

  _nominal_gait_period = _swing_time + _nominal_stance_time;
  _front_swing_time = _nominal_gait_period / 2. - 2. * dt;

  _front_previous_stance = _nominal_stance_time;
  _front_previous_swing = _swing_time;
  _front_current_stance = _nominal_stance_time;

  _hind_previous_stance = _nominal_stance_time;
  _hind_previous_swing = _swing_time;
  _hind_current_stance = _nominal_stance_time;

  _aerial_duration = 0.;

  _ini_fr = (WBCtrl::_model)._pGC[linkID::FR] -
            (WBCtrl::_model)._pGC[linkID::FR_abd];
  _ini_fl = (WBCtrl::_model)._pGC[linkID::FL] -
            (WBCtrl::_model)._pGC[linkID::FL_abd];

  _fin_fr = _ini_fr + _vel_des * _nominal_stance_time / 2.;
  _fin_fl = _ini_fl + _vel_des * _nominal_stance_time / 2.;

  _ini_hr = (WBCtrl::_model)._pGC[linkID::HR] -
            (WBCtrl::_model)._pGC[linkID::HR_abd];
  _ini_hl = (WBCtrl::_model)._pGC[linkID::HL] -
            (WBCtrl::_model)._pGC[linkID::HL_abd];

  _fin_hr = _ini_hr + _vel_des * _nominal_stance_time / 2.;
  _fin_hl = _ini_hl + _vel_des * _nominal_stance_time / 2.;

  _front_start_time = _curr_time;
  _hind_start_time = _curr_time;

  _ini_front_body = 0.5 * (WBCtrl::_model)._pGC[linkID::FR_abd] +
                    0.5 * (WBCtrl::_model)._pGC[linkID::FL_abd] -
                    0.5 * (WBCtrl::_model)._pGC[linkID::FR] -
                    0.5 * (WBCtrl::_model)._pGC[linkID::FL];

  _ini_hind_body = 0.5 * (WBCtrl::_model)._pGC[linkID::HR_abd] +
                   0.5 * (WBCtrl::_model)._pGC[linkID::HL_abd] -
                   0.5 * (WBCtrl::_model)._pGC[linkID::HR] -
                   0.5 * (WBCtrl::_model)._pGC[linkID::HL];
}

template <typename T>
void BoundingCtrl<T>::_contact_update() {
  typename std::vector<ContactSpec<T> *>::iterator iter =
      WBCtrl::_contact_list.begin();
  while (iter < WBCtrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
BoundingCtrl<T>::~BoundingCtrl() {  }

template <typename T>
void BoundingCtrl<T>::_lcm_data_sending() {}

template class BoundingCtrl<double>;
template class BoundingCtrl<float>;
