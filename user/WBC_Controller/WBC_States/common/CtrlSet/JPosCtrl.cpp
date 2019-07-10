#include "JPosCtrl.hpp"

#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/FixedBodyContact.hpp>

#include <Utilities/utilities.h>
#include <WBC/WBLC/WBLC.hpp>

template <typename T>
JPosCtrl<T>::JPosCtrl(const FloatingBaseModel<T>* robot)
    : Controller<T>(robot),
      _Kp(cheetah::num_act_joint),
      _Kd(cheetah::num_act_joint),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _des_jacc(cheetah::num_act_joint),
      _jpos_target(cheetah::num_act_joint),
      _end_time(1000.0),
      _dim_contact(0),
      _ctrl_start_time(0.) {
  _contact = new FixedBodyContact<T>();

  Ctrl::_contact_list.push_back(_contact);
  _dim_contact = Ctrl::_contact_list[0]->getDim();

  _wblc = new WBLC<T>(cheetah::dim_config, Ctrl::_contact_list);
  _wblc_data = new WBLC_ExtraData<T>();

  _wblc_data->W_qddot_ = DVec<T>::Constant(cheetah::dim_config, 100.0);
  _wblc_data->W_rf_ = DVec<T>::Constant(_dim_contact, 0.1);
  _wblc_data->W_xddot_ = DVec<T>::Constant(_dim_contact, 1000.0);

  // torque limit default setting
  _wblc_data->tau_min_ = DVec<T>::Constant(cheetah::num_act_joint, -150.);
  _wblc_data->tau_max_ = DVec<T>::Constant(cheetah::num_act_joint, 150.);

  _sp = StateProvider<T>::getStateProvider();

  printf("[Config JPos Control] Constructed\n");
}

template <typename T>
JPosCtrl<T>::~JPosCtrl() {
  delete _param_handler;
  delete _wblc;
  delete _wblc_data;

  typename std::vector<ContactSpec<T>*>::iterator iter2 =
      Ctrl::_contact_list.begin();

  while (iter2 < Ctrl::_contact_list.end()) {
    delete (*iter2);
    ++iter2;
  }
  Ctrl::_contact_list.clear();
}

template <typename T>
void JPosCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;

  DVec<T> gamma = DVec<T>::Zero(cheetah::num_act_joint);

  _contact_setup();
  _task_setup();
  _compute_torque_wblc(gamma);

  for (size_t leg(0); leg < cheetah::num_leg; ++leg) {
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
      ((LegControllerCommand<T>*)_cmd)[leg].tauFeedForward[jidx] =
          gamma[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qDes[jidx] =
          _des_jpos[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qdDes[jidx] =
          _des_jvel[cheetah::num_leg_joint * leg + jidx];
    }
  }
  Ctrl::_PostProcessing_Command();
}

template <typename T>
void JPosCtrl<T>::_compute_torque_wblc(DVec<T>& gamma) {
  // WBLC
  _wblc->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  DVec<T> des_jacc_cmd =
      _des_jacc + _Kp.cwiseProduct(_des_jpos - Ctrl::_robot_sys->_state.q) +
      _Kd.cwiseProduct(_des_jvel - Ctrl::_robot_sys->_state.qd);

  // pretty_print(des_jacc_cmd, std::cout, "des jacc");
  _wblc_data->_des_jacc_cmd = des_jacc_cmd;
  _wblc->MakeTorque(gamma, _wblc_data);
}

template <typename T>
void JPosCtrl<T>::_task_setup() {
  _des_jpos = _jpos_ini;
  _des_jvel.setZero();
  _des_jacc.setZero();

  if (_motion_type == MotionType::move_to_target) {
    for (size_t i(0); i < cheetah::num_act_joint; ++i) {
      _des_jpos[i] = smooth_change(_jpos_ini[i], _jpos_target[i], _end_time,
                                   Ctrl::_state_machine_time);
      _des_jvel[i] = smooth_change_vel(_jpos_ini[i], _jpos_target[i], _end_time,
                                       Ctrl::_state_machine_time);
      _des_jacc[i] = smooth_change_acc(_jpos_ini[i], _jpos_target[i], _end_time,
                                       Ctrl::_state_machine_time);
    }
  } else if (_motion_type == MotionType::swing) {
    static double vel_buffer(0.);
    if (Ctrl::_state_machine_time < 0.5) {
      vel_buffer = Ctrl::_state_machine_time / 0.5;
    } else {
      vel_buffer = 1.;
    }
    for (size_t i(0); i < cheetah::num_act_joint; ++i) {
      double omega(2. * M_PI * _swing_freq[i]);

      _des_jpos[i] = _swing_amp[i] * sin(omega * Ctrl::_state_machine_time);
      _des_jpos[i] += _jpos_target[i];

      _des_jvel[i] = vel_buffer * _swing_amp[i] * omega *
                     cos(omega * Ctrl::_state_machine_time);
      _des_jacc[i] = _swing_amp[i] * omega * omega *
                     sin(omega * Ctrl::_state_machine_time);
    }
  } else if (_motion_type == MotionType::stay) {
    _des_jpos = _jpos_target;
  }
}

template <typename T>
void JPosCtrl<T>::_contact_setup() {
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void JPosCtrl<T>::FirstVisit() {
  _jpos_ini = Ctrl::_robot_sys->_state.q;
  _ctrl_start_time = _sp->_curr_time;
}

template <typename T>
void JPosCtrl<T>::LastVisit() {}

template <typename T>
bool JPosCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > _end_time) {
    return true;
  }
  return false;
}

template <typename T>
void JPosCtrl<T>::CtrlInitialization(const std::string& category_name) {
  if (_param_handler) {
    std::string tmp_str;
    _param_handler->getString(category_name, "motion_type", tmp_str);
    if (tmp_str == "swing") {
      _motion_type = MotionType::swing;
      _param_handler->getVector<T>(category_name, "amplitude", _swing_amp);
      _param_handler->getVector<T>(category_name, "frequency", _swing_freq);
      _param_handler->getVector<T>(category_name, "phase", _swing_phase);
    } else if (tmp_str == "stay") {
      _motion_type = MotionType::stay;
    } else if (tmp_str == "move_to_target") {
      _motion_type = MotionType::move_to_target;
    } else {
      printf("[JPos Ctrl] Invalid Motion Type\n");
    }
    _param_handler->getValue<T>(category_name, "move_time", _end_time);
  } else {
    printf("[JPosCtrl] Yaml file is not properly set\n");
  }
}

template <typename T>
void JPosCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  _test_file_name = test_file;

  std::vector<T> tmp_vec;
  _param_handler->getVector<T>("target_jpos", tmp_vec);
  for (size_t i(0); i < cheetah::num_act_joint; ++i)
    _jpos_target[i] = tmp_vec[i];

  // Feedback Gain
  _param_handler->getVector<T>("Kp", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>("Kd", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _Kd[i] = tmp_vec[i];
  }
  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);
}

template class JPosCtrl<double>;
template class JPosCtrl<float>;
