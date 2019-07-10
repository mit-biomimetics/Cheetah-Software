#include "LegHeightTask.hpp"
// (X, Y, Z)
#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Utilities/Utilities_print.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
LegHeightTask<T>::LegHeightTask(const FloatingBaseModel<T>* robot, int link_idx,
                                int local_frame_idx)
    : Task<T>(1),
      _robot_sys(robot),
      _link_idx(link_idx),
      _local_frame_idx(local_frame_idx) {
  TK::Jt_ = DMat<T>::Zero(TK::dim_task_, cheetah::dim_config);
  TK::JtDotQdot_ = DVec<T>::Zero(TK::dim_task_);

  _Kp_kin = DVec<T>::Constant(TK::dim_task_, 1.);
  _Kp = DVec<T>::Constant(TK::dim_task_, 50.);
  _Kd = DVec<T>::Constant(TK::dim_task_, 3.);
}

template <typename T>
LegHeightTask<T>::~LegHeightTask() {}

template <typename T>
bool LegHeightTask<T>::_UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                                      const DVec<T>& acc_des) {
  T* pos_cmd = (T*)pos_des;
  Vec3<T> link_pos = _robot_sys->_pGC[_link_idx];

  Vec3<T> local_pos, local_vel;
  local_pos = _robot_sys->_pGC[_local_frame_idx];
  local_vel = _robot_sys->_vGC[_local_frame_idx];
  // Vec3<T> offset; offset.setZero();
  //_robot_sys->getPositionVelocity(_local_frame_idx, offset, local_pos,
  //local_vel);

  T leg_length = link_pos[2] - local_pos[2];
  // Z
  TK::pos_err_[0] = _Kp_kin[0] * ((*pos_cmd) - leg_length);
  TK::vel_des_[0] = vel_des[0];
  TK::acc_des_[0] = acc_des[0];

  TK::op_cmd_[0] = _Kp[0] * ((*pos_cmd) - leg_length) +
                   _Kd[0] * (TK::vel_des_[0] -
                             (_robot_sys->_vGC[_link_idx][2] - local_vel[2])) +
                   TK::acc_des_[0];

  // Quat<T> quat = _robot_sys->_state.bodyOrientation;
  // Mat3<T> Rot = ori::quaternionToRotationMatrix(quat);
  // TK::pos_err_ = Rot * TK::pos_err_;
  // TK::vel_des_ = Rot * TK::vel_des_;
  // TK::acc_des_ = Rot * TK::acc_des_;

  // printf("[Body Pos Task]\n");
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(TK::pos_err_, std::cout, "pos_err_");
  // pretty_print(*pos_cmd, std::cout, "pos cmd");
  // pretty_print(TK::Jt_, std::cout, "Jt");

  return true;
}

template <typename T>
bool LegHeightTask<T>::_UpdateTaskJacobian() {
  TK::Jt_ =
      (_robot_sys->_Jc[_link_idx]).block(2, 0, 1, cheetah::dim_config) -
      (_robot_sys->_Jc[_local_frame_idx]).block(2, 0, 1, cheetah::dim_config);
  // pretty_print(_robot_sys->_Jc[_link_idx], std::cout, "Jc link idx");
  // pretty_print(_robot_sys->_Jc[_local_frame_idx], std::cout, "Jc local");

  // TK::Jt_.block(0,3, 3,3) = Rot;
  // pretty_print(TK::Jt_, std::cout, "Jt");
  // TK::Jt_.block(0,3, 3,3) = Rot*TK::Jt_.block(0,3,3,3);
  return true;
}

template <typename T>
bool LegHeightTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class LegHeightTask<double>;
template class LegHeightTask<float>;
