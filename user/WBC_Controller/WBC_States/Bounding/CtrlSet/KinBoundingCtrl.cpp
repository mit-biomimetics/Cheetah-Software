#include "KinBoundingCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/Bounding/TaskSet/BodyRyRzTask.hpp>
#include <WBC_States/Bounding/TaskSet/BodyXYTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalHeadPosTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalPosTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalRollTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalTailPosTask.hpp>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>

#include <Utilities/save_file.h>
#include <ParamHandler/ParamHandler.hpp>
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/Bounding/BoundingTest.hpp>

template <typename T>
KinBoundingCtrl<T>::KinBoundingCtrl(BoundingTest<T>* bounding_test,
                                    const FloatingBaseModel<T>* robot)
    : Controller<T>(robot),
      _bounding_test(bounding_test),
      _jump_signal(false),
      _b_jump_initiation(false),
      _step_width(0.05),
      _contact_vel_threshold(3.0),
      _K_time(0.5),
      _front_foot_offset(0.02),
      _hind_foot_offset(-0.05),
      _swing_height(0.05),
      _end_time(10000.0),
      _dim_contact(0),
      _ctrl_start_time(0.),
      _step_length_lim(0.3),
      _fr_foot_vel(3),
      _fr_foot_acc(3),
      _fl_foot_vel(3),
      _fl_foot_acc(3),
      _hr_foot_vel(3),
      _hr_foot_acc(3),
      _hl_foot_vel(3),
      _hl_foot_acc(3),
      _total_mass(9.),
      _Fr_result(12),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _des_jacc(cheetah::num_act_joint),
      _wbcLCM(getLcmUrl(255)) {
  // Start from front swing & hind stance
  _b_front_swing = true;
  _b_hind_swing = false;

  //_b_front_swing = false;
  //_b_hind_swing = true;

  _fr_foot_vel.setZero();
  _fr_foot_acc.setZero();

  _fl_foot_vel.setZero();
  _fl_foot_acc.setZero();

  _hr_foot_vel.setZero();
  _hr_foot_acc.setZero();

  _hl_foot_vel.setZero();
  _hl_foot_acc.setZero();

  _local_roll_task = new LocalRollTask<T>(Ctrl::_robot_sys);
  _body_ryrz_task = new BodyRyRzTask<T>(Ctrl::_robot_sys);

  _fr_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::FR, linkID::FR_abd);
  _fl_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::FL, linkID::FL_abd);
  _hr_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::HR, linkID::HR_abd);
  _hl_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::HL, linkID::HL_abd);

  _local_head_pos_task = new LocalHeadPosTask<T>(Ctrl::_robot_sys);
  _local_tail_pos_task = new LocalTailPosTask<T>(Ctrl::_robot_sys);

  _fr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FR);
  _fl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FL);
  _hr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HR);
  _hl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HL);

  _kin_wbc = new KinWBC<T>(cheetah::dim_config);
  _wbic = new WBIC<T>(cheetah::dim_config, &(Ctrl::_contact_list),
                      &(Ctrl::_task_list));

  _wbic_data = new WBIC_ExtraData<T>();
  _wbic_data->_W_floating = DVec<T>::Constant(6, 5.);
  //_wbic_data->_W_floating = DVec<T>::Constant(6, 0.01);
  _sp = StateProvider<T>::getStateProvider();

  _folder_name = "/user/WBC_Controller/sim_data/";
  create_folder(_folder_name);

  _Fr_result.setZero();
  printf("[Kinematics Bounding Control] Constructed\n");
}

template <typename T>
void KinBoundingCtrl<T>::_ContactUpdate() {
  // Cheater Mode **************************************************************
  // T fr_z = Ctrl::_robot_sys->_pGC[linkID::FR][2];
  // T fl_z = Ctrl::_robot_sys->_pGC[linkID::FL][2];
  // T hr_z = Ctrl::_robot_sys->_pGC[linkID::HR][2];
  // T hl_z = Ctrl::_robot_sys->_pGC[linkID::HL][2];
  // T threshold(0.001);
  // if( (fr_z < threshold) && (fl_z < threshold) ){ _b_front_contact = true; }
  // else { _b_front_contact = false; }
  // if( (hr_z < threshold) && (hl_z < threshold) ){ _b_hind_contact = true; }
  // else { _b_hind_contact = false; }
  // END of Cheater Mode ******************************************************

  if (_b_front_swing) {
    _b_front_contact_est = false;
  } else {
    _b_front_contact_est = true;
  }

  if (_b_hind_swing) {
    _b_hind_contact_est = false;
  } else {
    _b_hind_contact_est = true;
  }

  DVec<T> qdot_pre = _qdot_pre_queue.front();
  T fr_vel_diff = _sp->_Qdot[6 + 2] - qdot_pre[6 + 2];
  T fl_vel_diff = _sp->_Qdot[6 + 5] - qdot_pre[6 + 5];
  T front_knee_vel_diff = std::max((fr_vel_diff), (fl_vel_diff));

  if (_b_front_swing && (_front_time > 0.6 * _swing_time)) {
    if (front_knee_vel_diff > _contact_vel_threshold ||
        (_front_time > 3.5 * _swing_time)) {
      _b_front_contact_est = true;

      // Jump initiation
      if(_jump_signal) {
        printf("jump initiation\n");
        _b_jump_initiation = true;
      }
    }
  }

  T hr_vel_diff = _sp->_Qdot[6 + 8] - qdot_pre[6 + 8];
  T hl_vel_diff = _sp->_Qdot[6 + 11] - qdot_pre[6 + 11];
  T hind_knee_vel_diff = std::max((hr_vel_diff), (hl_vel_diff));

  if (_b_hind_swing && (_hind_time > 0.6 * _swing_time)) {
    if ((hind_knee_vel_diff > _contact_vel_threshold) ||
        (_hind_time > 3.5 * _swing_time)) {
      _b_hind_contact_est = true;
    }
  }
}

template <typename T>
void KinBoundingCtrl<T>::_StatusCheck() {
  _stance_time = _nominal_stance_time;
  _gait_period = _nominal_gait_period;

  // High speed bounding stance time scaling
  if (fabs(_bounding_test->_body_vel[0]) * _nominal_stance_time >
      _step_length_lim) {
    _stance_time = fabs(_step_length_lim / _bounding_test->_body_vel[0]);
    _gait_period = _stance_time + _swing_time;

    // printf("vel cmd, step_length_lim, stance time: %f, %f, %f\n",
    //_bounding_test->_body_vel[0], _step_length_lim, _stance_time);
  }

  T adjust_time(0.);
  if (_b_front_swing && (_front_time > 0.5 * _swing_time)) {
    // Check Contact
    if (_b_front_contact_est) {
      _b_front_swing = false;
      _front_start_time = _sp->_curr_time;

      _front_previous_swing = _front_time;
      adjust_time = _K_time * (_hind_previous_stance + _aerial_duration -
                               0.5 * _gait_period);
      adjust_time = coerce<T>(adjust_time, 0, 0.04);
      _front_current_stance = _stance_time - adjust_time;

      T apex = _total_mass * 9.81 * (_swing_time + _front_current_stance) /
               (2. * 2.0 * 0.7 * _front_current_stance);
      apex *= _impact_amp;


      _front_z_impulse.setCurve(apex, _front_current_stance);
      _front_previous_stance = _front_current_stance;

      _ini_front_body = 0.5 * Ctrl::_robot_sys->_pGC[linkID::FR_abd] +
                        0.5 * Ctrl::_robot_sys->_pGC[linkID::FL_abd] -
                        0.5 * Ctrl::_robot_sys->_pGC[linkID::FR] -
                        0.5 * Ctrl::_robot_sys->_pGC[linkID::FL];

      _front_swing_time = _swing_time - 2. * Test<T>::dt;

      // Front time reset
      _front_time = 0.;
    }
  }

  if (_b_hind_swing && (_hind_time > 0.5 * _swing_time)) {
    // Check Contact
    if (_b_hind_contact_est) {
      _b_hind_swing = false;
      _hind_start_time = _sp->_curr_time;

      _hind_previous_swing = _hind_time;
      adjust_time = _K_time * (_front_previous_stance + _aerial_duration -
                               0.5 * _gait_period);
      adjust_time = coerce<T>(adjust_time, 0., 0.04);
      _hind_current_stance = _stance_time - adjust_time;

      T apex = _total_mass * 9.81 * (_swing_time + _hind_current_stance) /
               (2. * 2.0 * 0.7 * _hind_current_stance);
      apex *= _impact_amp;

      Vec3<T> rpy = ori::quatToRPY(Ctrl::_robot_sys->_state.bodyOrientation);
      apex *= (1. - _K_pitch * rpy[1]);

      _hind_z_impulse.setCurve(apex, _hind_current_stance);
      _hind_previous_stance = _hind_current_stance;

      _ini_hind_body = 0.5 * Ctrl::_robot_sys->_pGC[linkID::HR_abd] +
                       0.5 * Ctrl::_robot_sys->_pGC[linkID::HL_abd] -
                       0.5 * Ctrl::_robot_sys->_pGC[linkID::HR] -
                       0.5 * Ctrl::_robot_sys->_pGC[linkID::HL];

      // Hind time reset
      _hind_time = 0.;
    }
  }

  T scale(1.0);
  // If stance time is over switch to swing
  if ((!_b_front_swing) && (_front_time > (_front_current_stance))) {
    _b_front_swing = true;

    _ini_fr = Ctrl::_robot_sys->_pGC[linkID::FR] -
              Ctrl::_robot_sys->_pGC[linkID::FR_abd];
    _ini_fl = Ctrl::_robot_sys->_pGC[linkID::FL] -
              Ctrl::_robot_sys->_pGC[linkID::FL_abd];

    _fin_fr = _bounding_test->_body_vel * _stance_time / 2. * scale;
    _fin_fl = _bounding_test->_body_vel * _stance_time / 2. * scale;

    _fin_fr[0] += _front_foot_offset;
    _fin_fl[0] += _front_foot_offset;

    _front_start_time = _sp->_curr_time;
    _front_time = 0.;
  }
  // If stance time is over switch to swing
  if ((!_b_hind_swing) && (_hind_time > (_hind_current_stance))) {
    _b_hind_swing = true;

    _ini_hr = Ctrl::_robot_sys->_pGC[linkID::HR] -
              Ctrl::_robot_sys->_pGC[linkID::HR_abd];
    _ini_hl = Ctrl::_robot_sys->_pGC[linkID::HL] -
              Ctrl::_robot_sys->_pGC[linkID::HL_abd];

    _fin_hr = _bounding_test->_body_vel * _stance_time / 2. * scale;
    _fin_hl = _bounding_test->_body_vel * _stance_time / 2. * scale;

    _fin_hr[0] += _hind_foot_offset;
    _fin_hl[0] += _hind_foot_offset;

    _hind_start_time = _sp->_curr_time;
    _hind_time = 0.;
  }
}

template <typename T>
void KinBoundingCtrl<T>::_setupTaskAndContactList() {
  if (_b_front_swing) {
    _kin_task_list.push_back(_fl_foot_local_task);
    _kin_task_list.push_back(_fr_foot_local_task);

    Ctrl::_task_list.push_back(_fl_foot_local_task);
    Ctrl::_task_list.push_back(_fr_foot_local_task);
  } else {
    _kin_task_list.push_back(_local_head_pos_task);
    // Ctrl::_task_list.push_back(_local_head_pos_task);

    // Contact
    Ctrl::_contact_list.push_back(_fr_contact);
    Ctrl::_contact_list.push_back(_fl_contact);

    _kin_contact_list.push_back(_fr_contact);
    _kin_contact_list.push_back(_fl_contact);
  }

  if (_b_hind_swing) {
    _kin_task_list.push_back(_hr_foot_local_task);
    _kin_task_list.push_back(_hl_foot_local_task);

    Ctrl::_task_list.push_back(_hr_foot_local_task);
    Ctrl::_task_list.push_back(_hl_foot_local_task);
  } else {
    _kin_task_list.push_back(_local_tail_pos_task);
    // Ctrl::_task_list.push_back(_local_tail_pos_task);

    // Contact
    Ctrl::_contact_list.push_back(_hr_contact);
    Ctrl::_contact_list.push_back(_hl_contact);

    _kin_contact_list.push_back(_hr_contact);
    _kin_contact_list.push_back(_hl_contact);
  }
}

template <typename T>
void KinBoundingCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();

  
  // Initialize all
  Ctrl::_contact_list.clear();
  Ctrl::_task_list.clear();
  _kin_contact_list.clear();
  _kin_task_list.clear();

  // Update Time
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;

  if(Ctrl::_state_machine_time > 2.5){ 
    _jump_signal = true; 
  }
  _front_time = _sp->_curr_time - _front_start_time;
  _hind_time = _sp->_curr_time - _hind_start_time;

  // Update Contact
  _ContactUpdate();
  // Update Current Stance and Swing Status
  _StatusCheck();

  if ((!_b_front_swing) || (!_b_hind_swing)) {
    _aerial_duration = 0.;
    _kin_task_list.push_back(_local_roll_task);
    _kin_task_list.push_back(_body_ryrz_task);

    Ctrl::_task_list.push_back(_local_roll_task);
    Ctrl::_task_list.push_back(_body_ryrz_task);
  } else {
    _aerial_duration += Test<T>::dt;
  }

  _setupTaskAndContactList();

  DVec<T> gamma = DVec<T>::Zero(cheetah::num_act_joint);
  _contact_update();
  _body_task_setup();
  _leg_task_setup();
  _compute_torque_wbic(gamma);

  for (size_t leg(0); leg < cheetah::num_leg; ++leg) {
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
      ((LegControllerCommand<T>*)_cmd)[leg].tauFeedForward[jidx] =
          gamma[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qDes[jidx] =
          _des_jpos[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qdDes[jidx] =
          _des_jvel[cheetah::num_leg_joint * leg + jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kpJoint(jidx, jidx) =
          _Kp_joint[jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kdJoint(jidx, jidx) =
          _Kd_joint[jidx];
    }
  }

  _qdot_pre_queue.push(_sp->_Qdot);
  if (_qdot_pre_queue.size() > 4) _qdot_pre_queue.pop();
  Ctrl::_PostProcessing_Command();

  // LCM
  _lcm_data_sending();
  //_save_file();
}

template <typename T>
void KinBoundingCtrl<T>::_compute_torque_wbic(DVec<T>& gamma) {
  _kin_wbc->FindConfiguration(_sp->_Q, _kin_task_list, _kin_contact_list,
                              _des_jpos, _des_jvel, _des_jacc);

  // WBIC
  _wbic->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  _wbic->MakeTorque(gamma, _wbic_data);
}

template <typename T>
void KinBoundingCtrl<T>::_body_task_setup() {
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

  if (!_b_front_swing) {  // Stance
                          // for(size_t i(0); i<2; ++i){
    // pos_des[i] = _ini_front_body[i] +
    // _bounding_test->_body_vel[i]*_front_time; vel_des[i] =
    // _bounding_test->_body_vel[i];
    //}
    vel_des[0] = _bounding_test->_body_vel[0];
    vel_des[1] = _bounding_test->_body_vel[1];
    _local_head_pos_task->UpdateTask(&pos_des, vel_des, acc_des);
  }
  if (!_b_hind_swing) {  // Stance
                         // for(size_t i(0); i<2; ++i){
    // pos_des[i] = _ini_hind_body[i] + _bounding_test->_body_vel[i]*_hind_time;
    // vel_des[i] = _bounding_test->_body_vel[i];
    //}
    vel_des[0] = _bounding_test->_body_vel[0];
    vel_des[1] = _bounding_test->_body_vel[1];
    _local_tail_pos_task->UpdateTask(&pos_des, vel_des, acc_des);
  }
}

template <typename T>
void KinBoundingCtrl<T>::_leg_task_setup() {
  // for Z (height)
  T amp(_swing_height / 2.);
  T omega(2. * M_PI / _swing_time);
  T t;

  _fr_foot_pos.setZero();
  _fl_foot_pos.setZero();
  _hr_foot_pos.setZero();
  _hl_foot_pos.setZero();

  // T leg_height = -_target_leg_height;
  DVec<T> vel_des(3);
  vel_des.setZero();
  DVec<T> acc_des(3);
  acc_des.setZero();

  DVec<T> Fr_des(3);
  Fr_des.setZero();

  if (_b_front_swing) {
    t = _front_time;
    if (_front_time > _swing_time - 2. * Test<T>::dt) {
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

    _fr_foot_local_task->UpdateTask(&(_fr_foot_pos), _fr_foot_vel,
                                    _fr_foot_acc);
    _fl_foot_local_task->UpdateTask(&(_fl_foot_pos), _fl_foot_vel,
                                    _fl_foot_acc);
  } else {  // Front Leg Stance

    _wbic_data->_W_rf = DVec<T>::Constant(6, 1.0);

    _wbic_data->_W_rf[2] = 50.0;
    _wbic_data->_W_rf[5] = 50.0;

    Fr_des[2] = _front_z_impulse.getValue(_front_time);
    _fr_contact->setRFDesired(Fr_des);
    _fl_contact->setRFDesired(Fr_des);
  }

  if (_b_hind_swing) {
    t = _hind_time;
    if (_hind_time > _swing_time - 2. * Test<T>::dt) {
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
      _wbic_data->_W_rf = DVec<T>::Constant(6, 1.0);
      _wbic_data->_W_rf[2] = 50.0;
      _wbic_data->_W_rf[5] = 50.0;

      Fr_des[2] = _hind_z_impulse.getValue(_hind_time);
      _hr_contact->setRFDesired(Fr_des);
      _hl_contact->setRFDesired(Fr_des);
    } else {  // Front and Hind stance
      _wbic_data->_W_rf = DVec<T>::Constant(12, 1.0);
      _wbic_data->_W_rf[2] = 50.0;
      _wbic_data->_W_rf[5] = 50.0;
      _wbic_data->_W_rf[8] = 50.0;
      _wbic_data->_W_rf[11] = 50.0;

      Fr_des[2] = _front_z_impulse.getValue(_front_time);
      _fr_contact->setRFDesired(Fr_des);
      _fl_contact->setRFDesired(Fr_des);

      Fr_des[2] = _hind_z_impulse.getValue(_hind_time);
      _hr_contact->setRFDesired(Fr_des);
      _hl_contact->setRFDesired(Fr_des);
    }
  }
}

template class KinBoundingCtrl<double>;
template class KinBoundingCtrl<float>;
