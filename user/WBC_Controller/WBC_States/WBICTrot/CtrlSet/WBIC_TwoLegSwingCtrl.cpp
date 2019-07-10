#include "WBIC_TwoLegSwingCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <Utilities/utilities.h>
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/WBICTrot/WBICTrotTest.hpp>

template <typename T>
WBIC_TwoLegSwingCtrl<T>::WBIC_TwoLegSwingCtrl(WBICTrotTest<T>* test,
                                              const FloatingBaseModel<T>* robot,
                                              size_t cp1, size_t cp2)
    : Controller<T>(robot),
      _trot_test(test),
      _cp1(cp1),
      _cp2(cp2),
      _end_time(100.),
      _dim_contact(0),
      ctrl_start_time_(0.) {
  _body_pos_task = new BodyPosTask<T>(Ctrl::_robot_sys);
  _body_ori_task = new BodyOriTask<T>(Ctrl::_robot_sys);
  Ctrl::_task_list.push_back(_body_pos_task);
  Ctrl::_task_list.push_back(_body_ori_task);

  _cp_pos_task1 = new LinkPosTask<T>(Ctrl::_robot_sys, _cp1, false);
  _cp_pos_task2 = new LinkPosTask<T>(Ctrl::_robot_sys, _cp2, false);
  Ctrl::_task_list.push_back(_cp_pos_task1);
  Ctrl::_task_list.push_back(_cp_pos_task2);

  _foot_pos_ini1.setZero();
  _foot_pos_des1.setZero();
  _foot_vel_des1 = DVec<T>::Zero(3);
  _foot_acc_des1 = DVec<T>::Zero(3);

  _foot_pos_ini2.setZero();
  _foot_pos_des2.setZero();
  _foot_vel_des2 = DVec<T>::Zero(3);
  _foot_acc_des2 = DVec<T>::Zero(3);

  fr_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::FR);
  fl_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::FL);
  hr_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::HR);
  hl_contact_ = new SingleContact<T>(Ctrl::_robot_sys, linkID::HL);

  if (_cp1 == linkID::FR && _cp2 == linkID::HL) {
    Ctrl::_contact_list.push_back(fl_contact_);
    Ctrl::_contact_list.push_back(hr_contact_);
  } else {
    Ctrl::_contact_list.push_back(fr_contact_);
    Ctrl::_contact_list.push_back(hl_contact_);
  }

  for (size_t i(0); i < Ctrl::_contact_list.size(); ++i) {
    _dim_contact += Ctrl::_contact_list[i]->getDim();
  }

  DVec<T> Fr_des(3);
  Fr_des.setZero();
  Fr_des[2] = 9. * 9.81 / 2.;
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->setRFDesired(Fr_des);
    ++iter;
  }

  _kin_wbc = new KinWBC<T>(cheetah::dim_config);
  wbic_ = new WBIC<T>(cheetah::dim_config, &(Ctrl::_contact_list),
                      &(Ctrl::_task_list));
  wbic_data_ = new WBIC_ExtraData<T>();
  wbic_data_->_W_floating = DVec<T>::Constant(6, 1000.);
  wbic_data_->_W_rf = DVec<T>::Constant(_dim_contact, 1.);

  // if(_cp1 == linkID::HL || _cp2 == linkID::HL){
  // Ctrl::_contact_list.erase(Ctrl::_contact_list.begin() + 3);
  //}
  // if(_cp1 == linkID::HR || _cp2 == linkID::HR){
  // Ctrl::_contact_list.erase(Ctrl::_contact_list.begin() + 2);
  //}
  // if(_cp1 == linkID::FL || _cp2 == linkID::FL){
  // Ctrl::_contact_list.erase(Ctrl::_contact_list.begin() + 1);
  //}
  // if(_cp1 == linkID::FR || _cp2 == linkID::FR){
  // Ctrl::_contact_list.erase(Ctrl::_contact_list.begin());
  //}

  _sp = StateProvider<T>::getStateProvider();

  printf("[WBIC_Two Leg Swing Ctrl] Constructed\n");
}

template <typename T>
WBIC_TwoLegSwingCtrl<T>::~WBIC_TwoLegSwingCtrl() {
  delete wbic_;
  delete _kin_wbc;
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
void WBIC_TwoLegSwingCtrl<T>::OneStep(void* _cmd) {
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
void WBIC_TwoLegSwingCtrl<T>::_compute_torque_wbic(DVec<T>& gamma) {
  // WBIC
  wbic_->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  wbic_->MakeTorque(gamma, wbic_data_);
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::_task_setup() {
  _trot_test->_des_jpos.setZero();
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
  // Yaw control only
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

  // set Foot trajectory
  _GetSinusoidalSwingTrajectory(_foot_pos_ini1, _target_loc1,
                                Ctrl::_state_machine_time, _foot_pos_des1,
                                _foot_vel_des1, _foot_acc_des1);
  _GetSinusoidalSwingTrajectory(_foot_pos_ini2, _target_loc2,
                                Ctrl::_state_machine_time, _foot_pos_des2,
                                _foot_vel_des2, _foot_acc_des2);

  // Mat3<T> Rot = rpyToRotMat(_trot_test->_body_ori_rpy);
  //_foot_vel_des1.setZero(); _foot_acc_des1.setZero();
  // computeFootLoc(Rot.transpose(), _default_target_foot_loc_1, _step_time,
  // pos_des, vel_des,
  //_trot_test->_body_ang_vel, _foot_pos_des1);

  //_foot_vel_des2.setZero(); _foot_acc_des2.setZero();
  // computeFootLoc(Rot.transpose(), _default_target_foot_loc_2, _step_time,
  // pos_des, vel_des,
  //_trot_test->_body_ang_vel, _foot_pos_des2);

  // for Z (height)
  // T amp(_swing_height/2.);
  // T omega ( 2.*M_PI /_end_time );

  // T t = Ctrl::_state_machine_time;
  //_foot_pos_des1[2] = _foot_pos_ini1[2] + amp * (1-cos(omega * t));
  //_foot_vel_des1[2] = amp * omega * sin(omega * t);
  //_foot_acc_des1[2] = amp * omega * omega * cos(omega * t);

  //_foot_pos_des2[2] = _foot_pos_ini2[2] + amp * (1-cos(omega * t));
  //_foot_vel_des2[2] = amp * omega * sin(omega * t);
  //_foot_acc_des2[2] = amp * omega * omega * cos(omega * t);

  // Capture Point
  // if(false){
  if (true) {
    Vec3<T> global_body_vel;
    Vec3<T> local_body_vel;
    for (size_t i(0); i < 3; ++i) {
      local_body_vel[i] =
          Ctrl::_robot_sys->_state.bodyVelocity[i + 3];  // Local
    }
    Mat3<T> Rot_curr = ori::quaternionToRotationMatrix(
        Ctrl::_robot_sys->_state.bodyOrientation);
    global_body_vel = Rot_curr.transpose() * local_body_vel;

    for (size_t i(0); i < 2; ++i) {
      _foot_pos_des1[i] += sqrt(_target_body_height / 9.81) *
                           (global_body_vel[i] - _trot_test->_body_vel[i]);
      _foot_pos_des2[i] += sqrt(_target_body_height / 9.81) *
                           (global_body_vel[i] - _trot_test->_body_vel[i]);
    }
  }

  // Vec3<T> step_size = _foot_pos_des1 - _default_target_foot_loc_1;
  // for(size_t i(0); i<2; ++i){
  // if(step_size[i] > _step_size_max[i]){
  // step_size[i] = _step_size_max[i];
  //_foot_vel_des1[i] = 0.;
  //_foot_acc_des1[i] = 0.;
  //}
  // if(step_size[i] < _step_size_min[i]){
  // step_size[i] = _step_size_min[i];
  //_foot_vel_des1[i] = 0.;
  //_foot_acc_des1[i] = 0.;
  //}
  //}
  //_foot_pos_des1 = _default_target_foot_loc_1 + step_size;

  // step_size = _foot_pos_des2 - _default_target_foot_loc_2;
  // for(size_t i(0); i<2; ++i){
  // if(step_size[i] > _step_size_max[i]){
  // step_size[i] = _step_size_max[i];
  //_foot_vel_des2[i] = 0.;
  //_foot_acc_des2[i] = 0.;
  //}
  // if(step_size[i] < _step_size_min[i]){
  // step_size[i] = _step_size_min[i];
  //_foot_vel_des2[i] = 0.;
  //_foot_acc_des2[i] = 0.;
  //}
  //}
  //_foot_pos_des2 = _default_target_foot_loc_2 + step_size;
  // pretty_print(_step_size_min, std::cout, "min");
  // pretty_print(_step_size_max, std::cout, "max");
  // pretty_print(step_size, std::cout, "step size");
  // pretty_print(_foot_pos_des1, std::cout, "foot_pos_des1");

  _cp_pos_task1->UpdateTask(&(_foot_pos_des1), _foot_vel_des1, _foot_acc_des1);
  _cp_pos_task2->UpdateTask(&(_foot_pos_des2), _foot_vel_des2, _foot_acc_des2);

  _kin_wbc->FindConfiguration(_sp->_Q, Ctrl::_task_list, Ctrl::_contact_list,
                              _trot_test->_des_jpos, _trot_test->_des_jvel,
                              _trot_test->_des_jacc);

  // pretty_print(_trot_test->_des_jpos, std::cout, "des_jpos");
  // pretty_print(_trot_test->_des_jvel, std::cout, "des_jvel");
  // pretty_print(_trot_test->_des_jacc, std::cout, "des_jacc");
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::_GetSinusoidalSwingTrajectory(
    const Vec3<T>& ini, const Vec3<T>& fin, const T& t, Vec3<T>& pos_des,
    DVec<T>& vel_des, DVec<T>& acc_des) {
  for (size_t i(0); i < 2; ++i) {
    pos_des[i] = smooth_change(ini[i], fin[i], _end_time, t);
    vel_des[i] = smooth_change_vel(ini[i], fin[i], _end_time, t);
    acc_des[i] = smooth_change_acc(ini[i], fin[i], _end_time, t);
  }
  // for Z (height)
  T amp(_swing_height / 2.);
  T omega(2. * M_PI / _end_time);

  pos_des[2] = ini[2] + amp * (1 - cos(omega * t));
  vel_des[2] = amp * omega * sin(omega * t);
  acc_des[2] = amp * omega * omega * cos(omega * t);
  // TEST
  // acc_des[0] = 0.;
  // acc_des[1] = 0.;
  // acc_des[2] = 0.;
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::_contact_setup() {
  typename std::vector<ContactSpec<T>*>::iterator iter =
      Ctrl::_contact_list.begin();
  while (iter < Ctrl::_contact_list.end()) {
    (*iter)->UpdateContactSpec();
    ++iter;
  }
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::FirstVisit() {
  ctrl_start_time_ = _sp->_curr_time;
  _ini_body_pos = Ctrl::_robot_sys->_state.bodyPosition;

  _foot_pos_ini1 = Ctrl::_robot_sys->_pGC[_cp1];
  _foot_pos_ini2 = Ctrl::_robot_sys->_pGC[_cp2];

  _target_loc1 = _trot_test->_front_foot_loc;
  _target_loc2 = _trot_test->_hind_foot_loc;

  Vec3<T> cmd_vel = _trot_test->_body_vel;
  Vec3<T> next_body_pos = _ini_body_pos + cmd_vel * _end_time;

  Mat3<T> Rot = rpyToRotMat(_trot_test->_body_ori_rpy);
  computeFootLoc(Rot.transpose(), _default_target_foot_loc_1, _step_time,
                 next_body_pos, cmd_vel, _trot_test->_body_ang_vel,
                 _target_loc1);

  computeFootLoc(Rot.transpose(), _default_target_foot_loc_2, _step_time,
                 next_body_pos, cmd_vel, _trot_test->_body_ang_vel,
                 _target_loc2);

  _target_loc1 += _landing_offset;
  _target_loc2 += _landing_offset;

  //_SetBspline(_foot_pos_ini1, _target_loc1, _foot_traj_1);
  //_SetBspline(_foot_pos_ini2, _target_loc2, _foot_traj_2);

  // pretty_print(Rot, std::cout, "Rot");
  // pretty_print(_trot_test->_body_pos, std::cout, "commanded body_pos");
  // pretty_print(next_body_pos, std::cout, "nx body_pos");
  // pretty_print(_trot_test->_body_vel, std::cout, "body vel");
  // pretty_print(_trot_test->_body_ang_vel, std::cout, "body ang vel");
  // pretty_print(_target_loc1, std::cout, "target loc 1");
  // pretty_print(_target_loc2, std::cout, "target loc 2");
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::LastVisit() {
  _trot_test->_jpos_des_pre = _trot_test->_des_jpos;
}

template <typename T>
bool WBIC_TwoLegSwingCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > (_end_time - 2. * Test<T>::dt)) {
    return true;
  }
  return false;
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::CtrlInitialization(
    const std::string& category_name) {
  std::vector<T> tmp_vec;
  _param_handler->getVector<T>(category_name, "default_target_foot_location_1",
                               tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _default_target_foot_loc_1[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>(category_name, "default_target_foot_location_2",
                               tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _default_target_foot_loc_2[i] = tmp_vec[i];
  }
  // pretty_print(tmp_vec, "default target foot");
  _param_handler->getValue<T>(category_name, "swing_height", _swing_height);

  _param_handler->getVector<T>(category_name, "landing_offset", tmp_vec);
  // pretty_print(tmp_vec, "landing offset");
  for (size_t i(0); i < 3; ++i) _landing_offset[i] = tmp_vec[i];

  _param_handler->getVector<T>(category_name, "Kp_foot", tmp_vec);
  for (size_t i(0); i < 3; ++i) {
    ((LinkPosTask<T>*)_cp_pos_task1)->_Kp[i] = tmp_vec[i];
    ((LinkPosTask<T>*)_cp_pos_task2)->_Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>(category_name, "Kd_foot", tmp_vec);
  for (size_t i(0); i < 3; ++i) {
    ((LinkPosTask<T>*)_cp_pos_task1)->_Kd[i] = tmp_vec[i];
    ((LinkPosTask<T>*)_cp_pos_task2)->_Kd[i] = tmp_vec[i];
  }

  _param_handler->getVector<T>(category_name, "foot_step_min", tmp_vec);
  for (size_t i(0); i < 2; ++i) {
    _step_size_min[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>(category_name, "foot_step_max", tmp_vec);
  for (size_t i(0); i < 2; ++i) {
    _step_size_max[i] = tmp_vec[i];
  }
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);
  _param_handler->getValue<T>("body_height", _target_body_height);
  _param_handler->getValue<T>("swing_time", _end_time);

  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);

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

  _step_time = 0.;
  _step_time += _end_time;

  T tmp_value;
  _param_handler->getValue<T>("transition_time", tmp_value);
  _step_time += tmp_value;
  _step_time += tmp_value;

  _param_handler->getValue<T>("stance_time", tmp_value);
  _step_time += tmp_value;

  _param_handler->getValue<T>("step_time_ratio", tmp_value);
  _step_time *= tmp_value;
}

template <typename T>
void WBIC_TwoLegSwingCtrl<T>::computeFootLoc(
    const Mat3<T>& Rot, const Vec3<T>& shoulder, const T& step_time,
    const Vec3<T>& body_pos, const Vec3<T>& body_vel,
    const Vec3<T>& body_ang_vel, Vec3<T>& foot_loc) {
  foot_loc = body_pos + Rot * shoulder +
             step_time / 2. * (body_vel + body_ang_vel.cross(Rot * shoulder));
  foot_loc += _target_body_height / 9.81 * body_vel.cross(body_ang_vel);
  foot_loc[2] = 0.;
}

template class WBIC_TwoLegSwingCtrl<double>;
template class WBIC_TwoLegSwingCtrl<float>;
