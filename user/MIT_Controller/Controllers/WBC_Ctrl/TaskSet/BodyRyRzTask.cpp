#include "BodyRyRzTask.hpp"
// (Ry, Rz)

#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Dynamics/Quadruped.h>
#include <Math/orientation_tools.h>
#include <Utilities/Utilities_print.h>


template <typename T>
BodyRyRzTask<T>::BodyRyRzTask(const FloatingBaseModel<T>* robot)
    : Task<T>(2), _robot_sys(robot) {
  TK::Jt_ = DMat<T>::Zero(TK::dim_task_, cheetah::dim_config);
  TK::Jt_.block(0, 1, 2, 2).setIdentity();
  TK::JtDotQdot_ = DVec<T>::Zero(TK::dim_task_);

  _Kp_kin = DVec<T>::Constant(TK::dim_task_, 1.);
  _Kp = DVec<T>::Constant(TK::dim_task_, 50.);
  _Kd = DVec<T>::Constant(TK::dim_task_, 3.);
}

template <typename T>
BodyRyRzTask<T>::~BodyRyRzTask() {}

template <typename T>
bool BodyRyRzTask<T>::_UpdateCommand(const void* pos_des, const DVec<T>& vel_des,
                                     const DVec<T>& acc_des) {
  Quat<T>* ori_cmd = (Quat<T>*)pos_des;
  Quat<T> link_ori = (_robot_sys->_state.bodyOrientation);

  Quat<T> link_ori_inv;
  link_ori_inv[0] = link_ori[0];
  link_ori_inv[1] = -link_ori[1];
  link_ori_inv[2] = -link_ori[2];
  link_ori_inv[3] = -link_ori[3];
  // link_ori_inv /= link_ori.norm();

  // Quat<T> ori_err = ori::quatProduct(*ori_cmd, link_ori_inv);

  // implicit error definition
  Quat<T> ori_err = ori::quatProduct(link_ori_inv, *ori_cmd);
  if (ori_err[0] < 0.) {
    ori_err *= (-1.);
  }
  Vec3<T> ori_err_so3;
  ori::quaternionToso3(ori_err, ori_err_so3);
  SVec<T> curr_vel = _robot_sys->_state.bodyVelocity;

  // Rx, Ry, Rz
  for (int i(0); i < 2; ++i) {
    TK::pos_err_[i] = _Kp_kin[i] * ori_err_so3[i + 1];
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * ori_err_so3[i + 1] +
                     _Kd[i] * (TK::vel_des_[i] - curr_vel[i + 1]) +
                     TK::acc_des_[i];
  }
  // printf("[Body Ori Pitch Yaw Task]\n");
  // pretty_print(TK::pos_err_, std::cout, "pos_err_");
  // pretty_print(*ori_cmd, std::cout, "des_ori");
  // pretty_print(link_ori, std::cout, "curr_ori");
  // pretty_print(ori_err, std::cout, "quat_err");

  // pretty_print(link_ori_inv, std::cout, "ori_inv");
  // pretty_print(ori_err, std::cout, "ori_err");
  // pretty_print(*ori_cmd, std::cout, "cmd");
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(TK::Jt_, std::cout, "Jt");

  return true;
}

template <typename T>
bool BodyRyRzTask<T>::_UpdateTaskJacobian() {
  return true;
}

template <typename T>
bool BodyRyRzTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class BodyRyRzTask<double>;
template class BodyRyRzTask<float>;
