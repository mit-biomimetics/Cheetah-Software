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
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/JPosTask.hpp>

#include <Utilities/save_file.h>
#include <ParamHandler/ParamHandler.hpp>
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/Bounding/BoundingTest.hpp>

template <typename T>
void KinBoundingCtrl<T>::FirstVisit() {
  _qdot_pre_queue.push(_sp->_Qdot);
  T apex = _total_mass * 9.81 * (_swing_time + _nominal_stance_time) /
           (2. * 2.0 * 0.7 * _nominal_stance_time);

  _front_z_impulse.setCurve(apex, _nominal_stance_time);
  _hind_z_impulse.setCurve(apex, _nominal_stance_time);

  _nominal_gait_period = _swing_time + _nominal_stance_time;
  _front_swing_time = _nominal_gait_period / 2. - 2. * Test<T>::dt;

  _front_previous_stance = _nominal_stance_time;
  _front_previous_swing = _swing_time;
  _front_current_stance = _nominal_stance_time;

  _hind_previous_stance = _nominal_stance_time;
  _hind_previous_swing = _swing_time;
  _hind_current_stance = _nominal_stance_time;

  _aerial_duration = 0.;

  _ini_fr = Ctrl::_robot_sys->_pGC[linkID::FR] -
            Ctrl::_robot_sys->_pGC[linkID::FR_abd];
  _ini_fl = Ctrl::_robot_sys->_pGC[linkID::FL] -
            Ctrl::_robot_sys->_pGC[linkID::FL_abd];

  _fin_fr = _ini_fr + _bounding_test->_body_vel * _nominal_stance_time / 2.;
  _fin_fl = _ini_fl + _bounding_test->_body_vel * _nominal_stance_time / 2.;

  _ini_hr = Ctrl::_robot_sys->_pGC[linkID::HR] -
            Ctrl::_robot_sys->_pGC[linkID::HR_abd];
  _ini_hl = Ctrl::_robot_sys->_pGC[linkID::HL] -
            Ctrl::_robot_sys->_pGC[linkID::HL_abd];

  _fin_hr = _ini_hr + _bounding_test->_body_vel * _nominal_stance_time / 2.;
  _fin_hl = _ini_hl + _bounding_test->_body_vel * _nominal_stance_time / 2.;

  _front_start_time = _sp->_curr_time;
  _hind_start_time = _sp->_curr_time;

  _ctrl_start_time = _sp->_curr_time;

  _ini_front_body = 0.5 * Ctrl::_robot_sys->_pGC[linkID::FR_abd] +
                    0.5 * Ctrl::_robot_sys->_pGC[linkID::FL_abd] -
                    0.5 * Ctrl::_robot_sys->_pGC[linkID::FR] -
                    0.5 * Ctrl::_robot_sys->_pGC[linkID::FL];

  _ini_hind_body = 0.5 * Ctrl::_robot_sys->_pGC[linkID::HR_abd] +
                   0.5 * Ctrl::_robot_sys->_pGC[linkID::HL_abd] -
                   0.5 * Ctrl::_robot_sys->_pGC[linkID::HR] -
                   0.5 * Ctrl::_robot_sys->_pGC[linkID::HL];
}

template <typename T>
void KinBoundingCtrl<T>::_contact_update() {
  typename std::vector<ContactSpec<T> *>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }

  iter = _kin_contact_list.begin();
  while (iter < _kin_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

///////////////////////////////////////////////////////////////////////////////
//  PARAMETER setup
///////////////////////////////////////////////////////////////////////////////

template <typename T>
void KinBoundingCtrl<T>::CtrlInitialization(const std::string &category_name) {
  _param_handler->getValue<T>(category_name, "K_time", _K_time);
  _param_handler->getValue<T>(category_name, "K_pitch", _K_pitch);
  _param_handler->getValue<T>(category_name, "impact_amp", _impact_amp);

  _param_handler->getValue<T>(category_name, "front_foot_offset",
                              _front_foot_offset);
  _param_handler->getValue<T>(category_name, "hind_foot_offset",
                              _hind_foot_offset);

  _param_handler->getValue<T>(category_name, "swing_height", _swing_height);
  _param_handler->getValue<T>(category_name, "contact_vel_threshold",
                              _contact_vel_threshold);

  _param_handler->getValue<T>(category_name, "total_mass", _total_mass);
  _param_handler->getValue<T>(category_name, "step_length_limit",
                              _step_length_lim);
  _param_handler->getValue<T>(category_name, "step_width", _step_width);
}

template <typename T>
void KinBoundingCtrl<T>::SetTestParameter(const std::string &test_file) {
  _param_handler = new ParamHandler(test_file);
  std::vector<T> tmp_vec;
  if (!_param_handler->getValue<T>("leg_height", _target_leg_height)) {
    printf("[BoundingInitiate] No leg height\n");
  }

  _param_handler->getValue<T>("swing_time", _swing_time);
  _param_handler->getValue<T>("default_stance_time", _nominal_stance_time);
  _param_handler->getValue<T>("bounding_time", _end_time);

  _param_handler->getVector<T>("Kp_ryrz", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((BodyRyRzTask<T> *)_body_ryrz_task)->_Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>("Kd_ryrz", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((BodyRyRzTask<T> *)_body_ryrz_task)->_Kd[i] = tmp_vec[i];
  }

  // Head Pos
  _param_handler->getVector<T>("Kp_head", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((LocalHeadPosTask<T> *)_local_head_pos_task)->_Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>("Kd_head", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((LocalHeadPosTask<T> *)_local_head_pos_task)->_Kd[i] = tmp_vec[i];
  }

  // Tail Pos
  _param_handler->getVector<T>("Kp_tail", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((LocalTailPosTask<T> *)_local_tail_pos_task)->_Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>("Kd_tail", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    ((LocalTailPosTask<T> *)_local_tail_pos_task)->_Kd[i] = tmp_vec[i];
  }

  // Local Roll
  T tmp_value;
  _param_handler->getValue<T>("Kp_roll", tmp_value);
  ((LocalRollTask<T> *)_local_roll_task)->_Kp[0] = tmp_value;
  _param_handler->getValue<T>("Kd_roll", tmp_value);
  ((LocalRollTask<T> *)_local_roll_task)->_Kd[0] = tmp_value;

  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);
}

template <typename T>
void KinBoundingCtrl<T>::LastVisit() {}

template <typename T>
bool KinBoundingCtrl<T>::EndOfPhase() {
  if (_sp->_mode == 12) {
    _bounding_test->_stop = true;
    return true;
  }
  if(_jump_signal && _b_jump_initiation){
    _jump_signal = false;
    _b_jump_initiation = false;
    return true;
  }
  return false;
}

template <typename T>
KinBoundingCtrl<T>::~KinBoundingCtrl() {
  delete _kin_wbc;
  delete _wbic;
  delete _wbic_data;
  delete _param_handler;

  typename std::vector<Task<T> *>::iterator iter = Ctrl::_task_list.begin();
  while (iter < Ctrl::_task_list.end()) {
    delete (*iter);
    ++iter;
  }
  Ctrl::_task_list.clear();

  typename std::vector<ContactSpec<T> *>::iterator iter2 =
      Ctrl::_contact_list.begin();
  while (iter2 < Ctrl::_contact_list.end()) {
    delete (*iter2);
    ++iter2;
  }
  Ctrl::_contact_list.clear();
}

///////////////////////////////////////////////////////////////////////////////
//  DATA save
///////////////////////////////////////////////////////////////////////////////
template <typename T>
void KinBoundingCtrl<T>::_lcm_data_sending() {
  _wbc_data_lcm.contact_est[0] = _b_front_contact_est;
  _wbc_data_lcm.contact_est[1] = _b_hind_contact_est;

  for (size_t i(0); i < cheetah::num_act_joint; ++i) {
    _wbc_data_lcm.jpos_cmd[i] = _des_jpos[i];
    _wbc_data_lcm.jvel_cmd[i] = _des_jvel[i];
    _wbc_data_lcm.jacc_cmd[i] = _des_jacc[i];
  }

  _wbc_data_lcm.body_ori[3] = Ctrl::_robot_sys->_state.bodyOrientation[3];
  for (size_t i(0); i < 3; ++i) {
    _wbc_data_lcm.body_ori[i] = Ctrl::_robot_sys->_state.bodyOrientation[i];
    _wbc_data_lcm.body_ang_vel[i] = Ctrl::_robot_sys->_state.bodyVelocity[i];

    _wbc_data_lcm.fr_foot_pos_cmd[i] = _fr_foot_pos[i];
    _wbc_data_lcm.fl_foot_pos_cmd[i] = _fl_foot_pos[i];
    _wbc_data_lcm.hr_foot_pos_cmd[i] = _hr_foot_pos[i];
    _wbc_data_lcm.hl_foot_pos_cmd[i] = _hl_foot_pos[i];

    _wbc_data_lcm.fr_foot_vel_cmd[i] = _fr_foot_vel[i];
    _wbc_data_lcm.fl_foot_vel_cmd[i] = _fl_foot_vel[i];
    _wbc_data_lcm.hr_foot_vel_cmd[i] = _hr_foot_vel[i];
    _wbc_data_lcm.hl_foot_vel_cmd[i] = _hl_foot_vel[i];

    _wbc_data_lcm.fr_foot_local_pos[i] =
        _fr_foot_pos[i] - _fr_foot_local_task->getPosError()[i];
    _wbc_data_lcm.fl_foot_local_pos[i] =
        _fl_foot_pos[i] - _fl_foot_local_task->getPosError()[i];
    _wbc_data_lcm.hr_foot_local_pos[i] =
        _hr_foot_pos[i] - _hr_foot_local_task->getPosError()[i];
    _wbc_data_lcm.hl_foot_local_pos[i] =
        _hl_foot_pos[i] - _hl_foot_local_task->getPosError()[i];

    _wbc_data_lcm.fr_foot_local_vel[i] =
        (Ctrl::_robot_sys->_vGC[linkID::FR])[i] -
        (Ctrl::_robot_sys->_vGC[linkID::FR_abd])[i];
    _wbc_data_lcm.fl_foot_local_vel[i] =
        (Ctrl::_robot_sys->_vGC[linkID::FL])[i] -
        (Ctrl::_robot_sys->_vGC[linkID::FL_abd])[i];
    _wbc_data_lcm.hr_foot_local_vel[i] =
        (Ctrl::_robot_sys->_vGC[linkID::HR])[i] -
        (Ctrl::_robot_sys->_vGC[linkID::HR_abd])[i];
    _wbc_data_lcm.hl_foot_local_vel[i] =
        (Ctrl::_robot_sys->_vGC[linkID::HL])[i] -
        (Ctrl::_robot_sys->_vGC[linkID::HL_abd])[i];

    if ((!_b_front_swing) && (!_b_hind_swing)) {
      _wbc_data_lcm.fr_Fr_des[i] = _fr_contact->getRFDesired()[i];
      _wbc_data_lcm.fl_Fr_des[i] = _fl_contact->getRFDesired()[i];
      _wbc_data_lcm.hr_Fr_des[i] = _hr_contact->getRFDesired()[i];
      _wbc_data_lcm.hl_Fr_des[i] = _hl_contact->getRFDesired()[i];

      _wbc_data_lcm.fr_Fr[i] = _wbic_data->_Fr[i];
      _wbc_data_lcm.fl_Fr[i] = _wbic_data->_Fr[i + 3];
      _wbc_data_lcm.hr_Fr[i] = _wbic_data->_Fr[i + 6];
      _wbc_data_lcm.hl_Fr[i] = _wbic_data->_Fr[i + 9];
    } else if ((!_b_front_swing)) {
      _wbc_data_lcm.fr_Fr_des[i] = _fr_contact->getRFDesired()[i];
      _wbc_data_lcm.fl_Fr_des[i] = _fl_contact->getRFDesired()[i];

      _wbc_data_lcm.fr_Fr[i] = _wbic_data->_Fr[i];
      _wbc_data_lcm.fl_Fr[i] = _wbic_data->_Fr[i + 3];
    } else if ((!_b_hind_swing)) {
      _wbc_data_lcm.hr_Fr_des[i] = _hr_contact->getRFDesired()[i];
      _wbc_data_lcm.hl_Fr_des[i] = _hl_contact->getRFDesired()[i];

      _wbc_data_lcm.hr_Fr[i] = _wbic_data->_Fr[i];
      _wbc_data_lcm.hl_Fr[i] = _wbic_data->_Fr[i + 3];
    }
  }
  _wbcLCM.publish("wbc_lcm_data", &_wbc_data_lcm);
}

template <typename T>
void KinBoundingCtrl<T>::_save_file() {
  saveValue(Ctrl::_state_machine_time, _folder_name, "time");
  saveValue(_b_hind_contact, _folder_name, "hind_contact");
  saveValue(_b_front_contact, _folder_name, "front_contact");

  saveValue(_b_front_contact_est, _folder_name, "front_contact_est");
  saveValue(_b_hind_contact_est, _folder_name, "hind_contact_est");

  saveVector(_fr_foot_pos, _folder_name, "fr_foot_pos_cmd");
  saveVector(_fr_foot_vel, _folder_name, "fr_foot_vel_cmd");
  saveVector(_fr_foot_acc, _folder_name, "fr_foot_acc_cmd");

  saveVector(_fl_foot_pos, _folder_name, "fl_foot_pos_cmd");
  saveVector(_fl_foot_vel, _folder_name, "fl_foot_vel_cmd");
  saveVector(_fl_foot_acc, _folder_name, "fl_foot_acc_cmd");

  saveVector(_hr_foot_pos, _folder_name, "hr_foot_pos_cmd");
  saveVector(_hr_foot_vel, _folder_name, "hr_foot_vel_cmd");
  saveVector(_hr_foot_acc, _folder_name, "hr_foot_acc_cmd");

  saveVector(_hl_foot_pos, _folder_name, "hl_foot_pos_cmd");
  saveVector(_hl_foot_vel, _folder_name, "hl_foot_vel_cmd");
  saveVector(_hl_foot_acc, _folder_name, "hl_foot_acc_cmd");

  saveVector(_fr_foot_local_task->getPosError(), _folder_name, "fr_foot_err");
  saveVector(_fl_foot_local_task->getPosError(), _folder_name, "fl_foot_err");
  saveVector(_hr_foot_local_task->getPosError(), _folder_name, "hr_foot_err");
  saveVector(_hl_foot_local_task->getPosError(), _folder_name, "hl_foot_err");

  saveVector(_des_jpos, _folder_name, "jpos_cmd");
  saveVector(_des_jvel, _folder_name, "jvel_cmd");
  saveVector(_des_jacc, _folder_name, "jacc_cmd");

  saveVector(_sp->_Q, _folder_name, "config");
  saveVector(_sp->_Qdot, _folder_name, "qdot");

  _Fr_des = DVec<T>::Zero(12);
  _Fr_result = DVec<T>::Zero(12);

  if ((!_b_front_swing) && (!_b_hind_swing)) {  // Full Stance
    _Fr_des.head(3) = _fr_contact->getRFDesired();
    _Fr_des.segment(3, 3) = _fl_contact->getRFDesired();
    _Fr_des.segment(6, 3) = _hr_contact->getRFDesired();
    _Fr_des.segment(9, 3) = _hl_contact->getRFDesired();

    _Fr_result = _wbic_data->_Fr;
  } else if (!_b_front_swing) {  // Front Stance
    _Fr_des.head(3) = _fr_contact->getRFDesired();
    _Fr_des.segment(3, 3) = _fl_contact->getRFDesired();

    _Fr_result.head(6) = _wbic_data->_Fr;
  } else if (!_b_hind_swing) {  // Hind Stance
    _Fr_des.segment(6, 3) = _hr_contact->getRFDesired();
    _Fr_des.segment(9, 3) = _hl_contact->getRFDesired();

    _Fr_result.tail(6) = _wbic_data->_Fr;
  }
  saveVector(_Fr_des, _folder_name, "Fr_des");
  saveVector(_Fr_result, _folder_name, "Fr_result");
}

template class KinBoundingCtrl<double>;
template class KinBoundingCtrl<float>;
