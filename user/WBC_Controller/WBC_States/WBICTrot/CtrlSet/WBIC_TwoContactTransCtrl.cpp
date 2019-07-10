#include "WBIC_TwoContactTransCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/WBICTrot/WBICTrotTest.hpp>

template <typename T>
WBIC_TwoContactTransCtrl<T>::WBIC_TwoContactTransCtrl(
    WBICTrotTest<T>* test, const FloatingBaseModel<T>* robot, size_t cp1,
    size_t cp2, int transit_dir)
    : Controller<T>(robot),
      _trot_test(test),
      _cp1(cp1),
      _cp2(cp2),
      _transit_dir(transit_dir),
      b_set_height_target_(false),
      end_time_(100.),
      _dim_contact(0),
      ctrl_start_time_(0.) {
  _body_pos_task = new BodyPosTask<T>(Ctrl::_robot_sys);
  _body_ori_task = new BodyOriTask<T>(Ctrl::_robot_sys);

  Ctrl::_task_list.push_back(_body_pos_task);
  Ctrl::_task_list.push_back(_body_ori_task);

  // Pushback sequence is important !!
  fr_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::FR);
  fl_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::FL);
  hr_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::HR);
  hl_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::HL);

  Ctrl::_contact_list.push_back(fr_contact_);
  Ctrl::_contact_list.push_back(fl_contact_);
  Ctrl::_contact_list.push_back(hr_contact_);
  Ctrl::_contact_list.push_back(hl_contact_);

  DVec<T> Fr_des(3);
  Fr_des.setZero();
  Fr_des[2] = 9. * 9.81 / 4.;
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->setRFDesired(Fr_des);
    ++iter;
  }

  kin_wbc_ = new KinWBC<T>(cheetah::dim_config);
  wbic_ = new WBIC<T>(cheetah::dim_config, &(Ctrl::_contact_list),
                      &(Ctrl::_task_list));
  wbic_data_ = new WBIC_ExtraData<T>();

  for (size_t i(0); i < Ctrl::_contact_list.size(); ++i) {
    _dim_contact += Ctrl::_contact_list[i]->getDim();
  }

  wbic_data_->_W_floating = DVec<T>::Constant(6, 1000.);
  wbic_data_->_W_rf = DVec<T>::Constant(_dim_contact, 1.);

  _sp = StateProvider<T>::getStateProvider();

  printf("[WBIC_Two Contact Transition Ctrl] Constructed\n");
}

template <typename T>
WBIC_TwoContactTransCtrl<T>::~WBIC_TwoContactTransCtrl() {
  delete wbic_;
  delete kin_wbc_;
  delete wbic_data_;

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
void WBIC_TwoContactTransCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();
  Ctrl::_state_machine_time = _sp->_curr_time - ctrl_start_time_;

  DVec<T> gamma;
  _contact_setup();
  _task_setup();
  _compute_torque_wbic(gamma);

  for (size_t leg(0); leg < cheetah::num_leg; ++leg) {
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
      ((LegControllerCommand<T>*)_cmd)[leg].tauFeedForward[jidx] =
          gamma[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qDes[jidx] =
          _trot_test->_des_jpos[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qdDes[jidx] =
          _trot_test->_des_jvel[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].kpJoint(jidx, jidx) =
          _Kp_joint[jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kdJoint(jidx, jidx) =
          _Kd_joint[jidx];
    }
  }
  Ctrl::_PostProcessing_Command();
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::_compute_torque_wbic(DVec<T>& gamma) {
  // WBIC
  wbic_->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  wbic_->MakeTorque(gamma, wbic_data_);
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::_task_setup() {
  _trot_test->_des_jpos = Ctrl::_robot_sys->_state.q;
  _trot_test->_des_jvel.setZero();
  _trot_test->_des_jacc.setZero();

  // Calculate IK for a desired height and orientation.
  Vec3<T> pos_des;
  pos_des.setZero();
  DVec<T> vel_des(3);
  vel_des.setZero();
  DVec<T> acc_des(3);
  acc_des.setZero();
  Vec3<T> rpy_des;
  rpy_des.setZero();
  DVec<T> ang_vel_des(_body_ori_task->getDim());
  ang_vel_des.setZero();

  for (size_t i(0); i < 3; ++i) {
    pos_des[i] = _trot_test->_body_pos[i];
    vel_des[i] = _trot_test->_body_vel[i];
    acc_des[i] = _trot_test->_body_acc[i];
  }
  rpy_des[2] = _trot_test->_body_ori_rpy[2];
  ang_vel_des[2] = _trot_test->_body_ang_vel[2];

  _body_pos_task->UpdateTask(&(pos_des), vel_des, acc_des);

  // Set Desired Orientation
  Quat<T> des_quat;
  des_quat.setZero();
  des_quat = ori::rpyToQuat(rpy_des);

  DVec<T> ang_acc_des(_body_ori_task->getDim());
  ang_acc_des.setZero();
  _body_ori_task->UpdateTask(&(des_quat), ang_vel_des, ang_acc_des);

  kin_wbc_->FindConfiguration(_sp->_Q, Ctrl::_task_list, Ctrl::_contact_list,
                              _trot_test->_des_jpos, _trot_test->_des_jvel,
                              _trot_test->_des_jacc);

  if (_transit_dir < 0) {
    T alpha(Ctrl::_state_machine_time / end_time_);  // 0->1
    _trot_test->_des_jpos = alpha * _trot_test->_des_jpos +
                            (1. - alpha) * _trot_test->_jpos_des_pre;
  }
  // pretty_print(_trot_test->_des_jpos, std::cout, "_trot_test->_des_jpos);
  // pretty_print(_trot_test->_des_jvel, std::cout, "_trot_test->_des_jvel);
  // pretty_print(_trot_test->_des_jacc, std::cout, "_trot_test->_des_jacc);
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::_contact_setup() {
  T alpha =
      0.5 * (1 - cos(M_PI * Ctrl::_state_machine_time / end_time_));  // 0 -> 1
  T upper_lim(100.);

  if (_transit_dir >
      0) {  // Decrease reaction force & Increase full acceleration
    upper_lim = max_rf_z_ + alpha * (min_rf_z_ - max_rf_z_);
  } else {
    upper_lim = min_rf_z_ + alpha * (max_rf_z_ - min_rf_z_);
  }

  if (_cp1 == linkID::FR || _cp2 == linkID::FR) {
    _SetContact(0, upper_lim);
  }
  if (_cp1 == linkID::FL || _cp2 == linkID::FL) {
    _SetContact(1, upper_lim);
  }

  if (_cp1 == linkID::HR || _cp2 == linkID::HR) {
    _SetContact(2, upper_lim);
  }

  if (_cp1 == linkID::HL || _cp2 == linkID::HL) {
    _SetContact(3, upper_lim);
  }

  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::_SetContact(const size_t& cp_idx,
                                              const T& upper_lim) {
  ((SingleContact<T>*)Ctrl::_contact_list[cp_idx])->setMaxFz(upper_lim);
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::FirstVisit() {
  ctrl_start_time_ = _sp->_curr_time;
  ini_body_pos_ = Ctrl::_robot_sys->_state.bodyPosition;
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::LastVisit() {
  // printf("[ContactTransBody] End\n");
}

template <typename T>
bool WBIC_TwoContactTransCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > (end_time_ - 2. * Test<T>::dt)) {
    return true;
  }
  return false;
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::CtrlInitialization(
    const std::string& category_name) {
  // ParamHandler handler(CheetahConfigPath + setting_file_name + ".yaml");
  _param_handler->getValue<T>(category_name, "max_rf_z", max_rf_z_);
  _param_handler->getValue<T>(category_name, "min_rf_z", min_rf_z_);
}

template <typename T>
void WBIC_TwoContactTransCtrl<T>::SetTestParameter(
    const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  if (_param_handler->getValue<T>("body_height", _body_height_cmd)) {
    b_set_height_target_ = true;
  }
  _param_handler->getValue<T>("transition_time", end_time_);

  std::vector<T> tmp_vec;
  std::vector<T> tmp_Kp;
  std::vector<T> tmp_Kd;

  // Feedback gain for kinematic tasks
  _param_handler->getVector<T>("Kp_body_pos_kin", tmp_vec);
  _param_handler->getVector<T>("Kp_body_pos", tmp_Kp);
  _param_handler->getVector<T>("Kd_body_pos", tmp_Kd);
  for (size_t i(0); i < _body_pos_task->getDim(); ++i) {
    ((BodyPosTask<T>*)_body_pos_task)->_Kp_kin[i] = tmp_vec[i];
    ((BodyPosTask<T>*)_body_pos_task)->_Kp[i] = tmp_Kp[i];
    ((BodyPosTask<T>*)_body_pos_task)->_Kd[i] = tmp_Kd[i];
  }

  // Orientation (Kin, Kp, Kd)
  _param_handler->getVector<T>("Kp_body_ori_kin", tmp_vec);
  _param_handler->getVector<T>("Kp_body_ori", tmp_Kp);
  _param_handler->getVector<T>("Kd_body_ori", tmp_Kd);
  for (size_t i(0); i < _body_ori_task->getDim(); ++i) {
    ((BodyOriTask<T>*)_body_ori_task)->_Kp_kin[i] = tmp_vec[i];
    ((BodyOriTask<T>*)_body_ori_task)->_Kp[i] = tmp_Kp[i];
    ((BodyOriTask<T>*)_body_ori_task)->_Kd[i] = tmp_Kd[i];
  }

  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);
}

template class WBIC_TwoContactTransCtrl<double>;
template class WBIC_TwoContactTransCtrl<float>;
