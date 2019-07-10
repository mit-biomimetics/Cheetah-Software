#include "BodyXYTask.hpp"
// (X, Y)
#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Utilities/Utilities_print.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
BodyXYTask<T>::BodyXYTask(const FloatingBaseModel<T>* robot)
    : Task<T>(2), _robot_sys(robot) {
  TK::Jt_ = DMat<T>::Zero(TK::dim_task_, cheetah::dim_config);
  TK::Jt_(0, 3) = 1.;
  TK::Jt_(1, 4) = 1.;
  TK::JtDotQdot_ = DVec<T>::Zero(TK::dim_task_);

  _Kp_kin = DVec<T>::Constant(TK::dim_task_, 1.);

  _Kp = DVec<T>::Constant(TK::dim_task_, 50.);
  _Kd = DVec<T>::Constant(TK::dim_task_, 3.);
}

template <typename T>
BodyXYTask<T>::~BodyXYTask() {}

template <typename T>
bool BodyXYTask<T>::_UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                                   const DVec<T>& acc_des) {
  DVec<T>* pos_cmd = (DVec<T>*)pos_des;
  Vec3<T> link_pos = _robot_sys->_state.bodyPosition;

  Quat<T> quat = _robot_sys->_state.bodyOrientation;
  Mat3<T> Rot = ori::quaternionToRotationMatrix(quat);

  SVec<T> curr_vel = _robot_sys->_state.bodyVelocity;
  curr_vel.tail(3) = Rot.transpose() * curr_vel.tail(3);

  // X, Y, Z
  for (size_t i(0); i < TK::dim_task_; ++i) {
    TK::pos_err_[i] = _Kp_kin[i] * ((*pos_cmd)[i] - link_pos[i]);
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * ((*pos_cmd)[i] - link_pos[i]) +
                     _Kd[i] * (TK::vel_des_[i] - curr_vel[i + 3]) +
                     TK::acc_des_[i];
  }

  // printf("[Body Pos Task]\n");
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(TK::pos_err_, std::cout, "pos_err_");
  // pretty_print(*pos_cmd, std::cout, "pos cmd");
  // pretty_print(TK::Jt_, std::cout, "Jt");

  return true;
}

template <typename T>
bool BodyXYTask<T>::_UpdateTaskJacobian() {
  Quat<T> quat = _robot_sys->_state.bodyOrientation;
  Mat3<T> Rot = ori::quaternionToRotationMatrix(quat);
  TK::Jt_.block(0, 3, 2, 3) = (Rot.transpose()).block(0, 0, 2, 3);
  // TK::Jt_.block(0,3, 3,3) = Rot;
  // pretty_print(TK::Jt_, std::cout, "Jt");
  // TK::Jt_.block(0,3, 3,3) = Rot*TK::Jt_.block(0,3,3,3);
  return true;
}

template <typename T>
bool BodyXYTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class BodyXYTask<double>;
template class BodyXYTask<float>;
