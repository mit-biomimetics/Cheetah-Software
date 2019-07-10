#include "FullContactTransCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <WBC/WBLC/KinWBC.hpp>
#include <WBC/WBLC/WBLC.hpp>

template <typename T>
FullContactTransCtrl<T>::FullContactTransCtrl(const FloatingBaseModel<T>* robot)
    : Controller<T>(robot),
      _Kp(cheetah::num_act_joint),
      _Kd(cheetah::num_act_joint),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _des_jacc(cheetah::num_act_joint),
      _b_set_height_target(false),
      _end_time(100.),
      _dim_contact(0),
      _ctrl_start_time(0.) {
  _body_pos_task = new BodyPosTask<T>(Ctrl::_robot_sys);
  _body_ori_task = new BodyOriTask<T>(Ctrl::_robot_sys);

  Ctrl::_task_list.push_back(_body_ori_task);
  Ctrl::_task_list.push_back(_body_pos_task);

  _fr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FR);
  _fl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FL);
  _hr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HR);
  _hl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HL);

  Ctrl::_contact_list.push_back(_fr_contact);
  Ctrl::_contact_list.push_back(_fl_contact);
  Ctrl::_contact_list.push_back(_hr_contact);
  Ctrl::_contact_list.push_back(_hl_contact);

  _kin_wbc = new KinWBC<T>(cheetah::dim_config);
  _wblc = new WBLC<T>(cheetah::dim_config, Ctrl::_contact_list);
  _wblc_data = new WBLC_ExtraData<T>();

  for (size_t i(0); i < Ctrl::_contact_list.size(); ++i) {
    _dim_contact += Ctrl::_contact_list[i]->getDim();
  }

  _wblc_data->W_qddot_ = DVec<T>::Constant(cheetah::dim_config, 100.0);
  _wblc_data->W_rf_ = DVec<T>::Constant(_dim_contact, 1.);
  _wblc_data->W_xddot_ = DVec<T>::Constant(_dim_contact, 1000.0);

  int idx_offset(0);
  for (size_t i(0); i < Ctrl::_contact_list.size(); ++i) {
    _wblc_data->W_rf_[idx_offset + Ctrl::_contact_list[i]->getFzIndex()] = 0.01;
    idx_offset += Ctrl::_contact_list[i]->getDim();
  }

  // torque limit default setting
  _wblc_data->tau_min_ = DVec<T>::Constant(cheetah::num_act_joint, -150.);
  _wblc_data->tau_max_ = DVec<T>::Constant(cheetah::num_act_joint, 150.);

  _sp = StateProvider<T>::getStateProvider();

  printf("[Full Contact Transition Ctrl] Constructed\n");
}

template <typename T>
FullContactTransCtrl<T>::~FullContactTransCtrl() {
  delete _wblc;
  delete _kin_wbc;
  delete _wblc_data;

  typename std::vector<Task<T>*>::iterator iter = Ctrl::_task_list.begin();
  while (iter < Ctrl::_task_list.end()) {
    delete (*iter);
    ++iter;
  }
  Ctrl::_task_list.clear();

  typename std::vector<ContactSpec<T>*>::iterator iter2 =
      Ctrl::_contact_list.begin();
  while (iter2 < Ctrl::_contact_list.end()) {
    delete (*iter2);
    ++iter2;
  }
  Ctrl::_contact_list.clear();
}

template <typename T>
void FullContactTransCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;

  DVec<T> gamma;
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

      ((LegControllerCommand<T>*)_cmd)[leg].kpJoint(jidx, jidx) =
          _Kp_joint[jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kdJoint(jidx, jidx) =
          _Kd_joint[jidx];
    }
  }
  Ctrl::_PostProcessing_Command();
}

template <typename T>
void FullContactTransCtrl<T>::_compute_torque_wblc(DVec<T>& gamma) {
  // WBLC
  _wblc->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  _wblc_data->_des_jacc_cmd =
      _des_jacc + _Kp.cwiseProduct(_des_jpos - Ctrl::_robot_sys->_state.q) +
      _Kd.cwiseProduct(_des_jvel - Ctrl::_robot_sys->_state.qd);

  _wblc->MakeTorque(gamma, _wblc_data);
  // pretty_print(_wblc_data->Fr_, std::cout, "fr full contact");
}

template <typename T>
void FullContactTransCtrl<T>::_task_setup() {
  _des_jpos = _ini_jpos;
  _des_jvel.setZero();
  _des_jacc.setZero();

  // Calculate IK for a desired height and orientation.
  if (!_b_set_height_target) {
    printf("No Height Command\n");
    exit(0);
  }
  Vec3<T> des_pos;
  des_pos.setZero();
  DVec<T> vel_des(3);
  vel_des.setZero();
  DVec<T> acc_des(3);
  acc_des.setZero();

  des_pos[2] = _ini_body_pos[2] + Ctrl::_state_machine_time / _end_time *
                                      (_target_body_height - _ini_body_pos[2]);

  _body_pos_task->UpdateTask(&(des_pos), vel_des, acc_des);

  // Set Desired Orientation
  Quat<T> des_quat;
  des_quat.setZero();
  des_quat[0] = 1.;

  DVec<T> ang_vel_des(_body_ori_task->getDim());
  ang_vel_des.setZero();
  DVec<T> ang_acc_des(_body_ori_task->getDim());
  ang_acc_des.setZero();
  _body_ori_task->UpdateTask(&(des_quat), ang_vel_des, ang_acc_des);

  _kin_wbc->FindConfiguration(_sp->_Q, Ctrl::_task_list, Ctrl::_contact_list,
                              _des_jpos, _des_jvel, _des_jacc);

  // pretty_print(_des_jpos, std::cout, "des_jpos");
  // pretty_print(_des_jvel, std::cout, "des_jvel");
  // pretty_print(_des_jacc, std::cout, "des_jacc");
}

template <typename T>
void FullContactTransCtrl<T>::_contact_setup() {
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    ((SingleContact<T>*)(*iter))
        ->setMaxFz(_min_rf_z + Ctrl::_state_machine_time / _end_time *
                                   (_max_rf_z - _min_rf_z));
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void FullContactTransCtrl<T>::FirstVisit() {
  _ini_jpos = Ctrl::_robot_sys->_state.q;
  _ctrl_start_time = _sp->_curr_time;
  _ini_body_pos = Ctrl::_robot_sys->_state.bodyPosition;
}

template <typename T>
void FullContactTransCtrl<T>::LastVisit() {
  // printf("[ContactTransBody] End\n");
}

template <typename T>
bool FullContactTransCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > _end_time) {
    return true;
  }
  return false;
}

template <typename T>
void FullContactTransCtrl<T>::CtrlInitialization(
    const std::string& category_name) {
  _param_handler->getValue<T>(category_name, "max_rf_z", _max_rf_z);
  _param_handler->getValue<T>(category_name, "min_rf_z", _min_rf_z);
  _param_handler->getValue<T>(category_name, "body_lifting_time", _end_time);
}

template <typename T>
void FullContactTransCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  if (_param_handler->getValue<T>("body_height", _target_body_height)) {
    _b_set_height_target = true;
  }

  // Feedback Gain
  std::vector<T> tmp_vec;
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
template class FullContactTransCtrl<double>;
template class FullContactTransCtrl<float>;
