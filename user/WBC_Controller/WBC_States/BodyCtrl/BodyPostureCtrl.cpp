#include "BodyPostureCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <ParamHandler/ParamHandler.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC/WBLC/WBLC.hpp>
#include <WBC_States/BodyCtrl/BodyCtrlTest.hpp>

template <typename T>
BodyPostureCtrl<T>::BodyPostureCtrl(const FloatingBaseModel<T>* robot,
                                    BodyCtrlTest<T>* test)
    : Controller<T>(robot),
      _Kp(cheetah::num_act_joint),
      _Kd(cheetah::num_act_joint),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _des_jacc(cheetah::num_act_joint),
      _end_time(1000.0),
      _dim_contact(0),
      _ctrl_start_time(0.) {
  _body_pos_task = new BodyPosTask<T>(Ctrl::_robot_sys);
  _body_ori_task = new BodyOriTask<T>(Ctrl::_robot_sys);

  Ctrl::_task_list.push_back(_body_pos_task);
  Ctrl::_task_list.push_back(_body_ori_task);

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
  _target_ori_command.setZero();

  _body_test = test;
  for (size_t i(0); i < 3; ++i) {
    _ori_cmd_filter.push_back(
        new digital_lp_filter<T>(2. * M_PI * 10., _body_test->dt));
  }
  printf("[Body Control] Constructed\n");
}

template <typename T>
BodyPostureCtrl<T>::~BodyPostureCtrl() {
  delete _wblc;
  delete _wblc_data;
  delete _kin_wbc;
  delete _param_handler;

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
void BodyPostureCtrl<T>::OneStep(void* _cmd) {
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

      ((LegControllerCommand<T>*)_cmd)[leg].kpJoint(jidx, jidx) =
          _Kp_joint[jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kdJoint(jidx, jidx) =
          _Kd_joint[jidx];
    }
  }
  Ctrl::_PostProcessing_Command();
}

template <typename T>
void BodyPostureCtrl<T>::_compute_torque_wblc(DVec<T>& gamma) {
  // WBLC
  _wblc->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  DVec<T> des_jacc_cmd =
      _des_jacc + _Kp.cwiseProduct(_des_jpos - Ctrl::_robot_sys->_state.q) +
      _Kd.cwiseProduct(_des_jvel - Ctrl::_robot_sys->_state.qd);

  _wblc_data->_des_jacc_cmd = des_jacc_cmd;
  _wblc->MakeTorque(gamma, _wblc_data);

  // pretty_print(Ctrl::_grav, std::cout, "grav");
  // pretty_print(Ctrl::_coriolis, std::cout, "coriolis");
  // pretty_print(Ctrl::_A, std::cout, "A");
}

template <typename T>
void BodyPostureCtrl<T>::_task_setup() {
  _des_jpos.setZero();
  _des_jvel.setZero();
  _des_jacc.setZero();

  // Calculate IK for a desired height and orientation.
  Vec3<T> pos_des;
  pos_des.setZero();
  DVec<T> vel_des(3);
  vel_des.setZero();
  DVec<T> acc_des(3);
  acc_des.setZero();
  pos_des[2] = _target_body_height;

  _body_pos_task->UpdateTask(&(pos_des), vel_des, acc_des);

  // Set Desired Orientation
  Quat<T> des_quat;
  des_quat.setZero();

  for (size_t i(0); i < 3; ++i) {
    //_ori_cmd_filter[i]->input(_sp->_ori_command[i]);
    //_target_ori_command[i] += _ori_cmd_filter[i]->output()*Test<T>::dt;

    _target_ori_command[i] = _sp->_ori_command[i];
  }

  // static int count(0);
  // count++;
  // if(count%500 ==0){
  // pretty_print(_target_ori_command, std::cout, "target ori");
  //}
  // Mat3<T> Rot = rpyToRotMat(_target_ori_command);
  // Eigen::Quaternion<T> eigen_quat(Rot.transpose());
  // des_quat[0] = eigen_quat.w();
  // des_quat[1] = eigen_quat.x();
  // des_quat[2] = eigen_quat.y();
  // des_quat[3] = eigen_quat.z();

  des_quat = rpyToQuat(_target_ori_command);

  DVec<T> ang_vel_des(_body_ori_task->getDim());
  ang_vel_des.setZero();
  DVec<T> ang_acc_des(_body_ori_task->getDim());
  ang_acc_des.setZero();
  _body_ori_task->UpdateTask(&(des_quat), ang_vel_des, ang_acc_des);

  _kin_wbc->FindConfiguration(_sp->_Q, Ctrl::_task_list, Ctrl::_contact_list,
                              _des_jpos, _des_jvel, _des_jacc);
}

template <typename T>
void BodyPostureCtrl<T>::_contact_setup() {
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void BodyPostureCtrl<T>::FirstVisit() {
  _ctrl_start_time = _sp->_curr_time;
}

template <typename T>
void BodyPostureCtrl<T>::LastVisit() {}

template <typename T>
bool BodyPostureCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > _end_time) {
    return true;
  }
  return false;
}

template <typename T>
void BodyPostureCtrl<T>::CtrlInitialization(const std::string& category_name) {
  (void)category_name;
}

template <typename T>
void BodyPostureCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  if (!_param_handler->getValue<T>("body_height", _target_body_height)) {
    printf("[Body Posture Ctrl] Error. No Height Target\n");
  }
  _param_handler->getValue<T>("stance_time", _end_time);

  std::vector<T> tmp_vec;
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

template class BodyPostureCtrl<double>;
template class BodyPostureCtrl<float>;
