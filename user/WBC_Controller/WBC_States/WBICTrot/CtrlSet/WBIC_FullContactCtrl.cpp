#include "WBIC_FullContactCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <ParamHandler/ParamHandler.hpp>
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/WBICTrot/WBICTrotTest.hpp>

template <typename T>
WBIC_FullContactCtrl<T>::WBIC_FullContactCtrl(WBICTrotTest<T>* trot_test,
                                              const FloatingBaseModel<T>* robot)
    : Controller<T>(robot),
      _trot_test(trot_test),
      _Kp(cheetah::num_act_joint),
      _Kd(cheetah::num_act_joint),
      _end_time(1000.0),
      _dim_contact(0),
      _ctrl_start_time(0.) {
  _body_ori_task = new BodyOriTask<T>(Ctrl::_robot_sys);
  _body_pos_task = new BodyPosTask<T>(Ctrl::_robot_sys);

  Ctrl::_task_list.push_back(_body_pos_task);
  Ctrl::_task_list.push_back(_body_ori_task);

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

  printf("[Body Control] Constructed\n");
}

template <typename T>
void WBIC_FullContactCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;

  DVec<T> gamma = DVec<T>::Zero(cheetah::num_act_joint);
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

  // static int count(0);
  //++count;
  // if(count>10){exit(0); }
}

template <typename T>
void WBIC_FullContactCtrl<T>::_compute_torque_wbic(DVec<T>& gamma) {
  // WBIC
  wbic_->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  wbic_->MakeTorque(gamma, wbic_data_);
}

template <typename T>
void WBIC_FullContactCtrl<T>::_task_setup() {
  _trot_test->_des_jpos.setZero();
  _trot_test->_des_jvel.setZero();
  _trot_test->_des_jacc.setZero();

  // Calculate IK for a desired height and orientation.
  Vec3<T> pos_des;
  pos_des.setZero();
  pos_des[2] = _trot_test->_body_pos[2];
  DVec<T> vel_des(3);
  vel_des.setZero();
  DVec<T> acc_des(3);
  acc_des.setZero();
  Vec3<T> rpy_des;
  rpy_des.setZero();
  rpy_des[2] = _ini_rpy[2];
  DVec<T> ang_vel_des(_body_ori_task->getDim());
  ang_vel_des.setZero();

  if (_sp->_mode == 11) {
    for (size_t i(0); i < 3; ++i) {
      pos_des[i] = _trot_test->_body_pos[i];
      vel_des[i] = _trot_test->_body_vel[i];
      acc_des[i] = _trot_test->_body_acc[i];
    }
    // Yaw control only
    rpy_des[2] = _trot_test->_body_ori_rpy[2];
    ang_vel_des[2] = _trot_test->_body_ang_vel[2];
  } else {
    _trot_test->_body_ori_rpy[2] = _ini_rpy[2];
  }
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
}

template <typename T>
void WBIC_FullContactCtrl<T>::_contact_setup() {
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void WBIC_FullContactCtrl<T>::FirstVisit() {
  _ctrl_start_time = _sp->_curr_time;
  _ini_rpy = _trot_test->_body_ori_rpy;
}

template <typename T>
void WBIC_FullContactCtrl<T>::LastVisit() {}

template <typename T>
bool WBIC_FullContactCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > (_end_time - 2. * Test<T>::dt) &&
      (_sp->_mode == 11)) {
    return true;
  }
  return false;
}

template <typename T>
void WBIC_FullContactCtrl<T>::CtrlInitialization(
    const std::string& category_name) {
  (void)category_name;
}

template <typename T>
void WBIC_FullContactCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  _param_handler->getValue<T>("stance_time", _end_time);

  std::vector<T> tmp_vec;
  std::vector<T> tmp_Kp;
  std::vector<T> tmp_Kd;

  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);

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
}

template <typename T>
WBIC_FullContactCtrl<T>::~WBIC_FullContactCtrl() {
  delete wbic_;
  delete wbic_data_;
  delete kin_wbc_;
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

template class WBIC_FullContactCtrl<double>;
template class WBIC_FullContactCtrl<float>;
