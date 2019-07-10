#include "LocalTailPosTask.hpp"

#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Dynamics/Quadruped.h>
#include <Math/orientation_tools.h>
#include <Utilities/Utilities_print.h>

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
LocalTailPosTask<T>::LocalTailPosTask(const FloatingBaseModel<T>* robot)
    : Task<T>(3), _robot_sys(robot) {
  TK::Jt_ = DMat<T>::Zero(TK::dim_task_, cheetah::dim_config);
  TK::JtDotQdot_ = DVec<T>::Zero(TK::dim_task_);

  _Kp_kin = DVec<T>::Constant(TK::dim_task_, 1.);
  _Kp = DVec<T>::Constant(TK::dim_task_, 50.);
  _Kd = DVec<T>::Constant(TK::dim_task_, 3.);
}

template <typename T>
LocalTailPosTask<T>::~LocalTailPosTask() {}

template <typename T>
bool LocalTailPosTask<T>::_UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                                         const DVec<T>& acc_des) {
  Vec3<T>* pos_cmd = (Vec3<T>*)pos_des;
  Vec3<T> link_pos, link_vel, local_pos, local_vel;
  link_pos = 0.5 * _robot_sys->_pGC[linkID::HR_abd] +
             0.5 * _robot_sys->_pGC[linkID::HL_abd];
  link_vel = 0.5 * _robot_sys->_vGC[linkID::HR_abd] +
             0.5 * _robot_sys->_vGC[linkID::HL_abd];

  local_pos =
      0.5 * _robot_sys->_pGC[linkID::HR] + 0.5 * _robot_sys->_pGC[linkID::HL];
  local_vel =
      0.5 * _robot_sys->_vGC[linkID::HR] + 0.5 * _robot_sys->_vGC[linkID::HL];

  // X, Y, Z
  for (size_t i(0); i < TK::dim_task_; ++i) {
    TK::pos_err_[i] =
        _Kp_kin[i] * ((*pos_cmd)[i] - (link_pos[i] - local_pos[i]));
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * ((*pos_cmd)[i] - (link_pos[i] - local_pos[i])) +
                     _Kd[i] * (link_vel[i] - local_vel[i]) + TK::acc_des_[i];
  }
  // printf("[Tail Local Position Task]\n");
  // pretty_print(TK::pos_err_, std::cout, "pos_err");
  // pretty_print(*pos_cmd, std::cout, "pos cmd");
  // pretty_print(link_pos, std::cout, "link pos");
  // pretty_print(local_pos, std::cout, "local pos");
  // pretty_print(ori_err, std::cout, "quat_err");

  // pretty_print(link_ori_inv, std::cout, "ori_inv");
  // pretty_print(ori_err, std::cout, "ori_err");
  // pretty_print(*ori_cmd, std::cout, "cmd");
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(TK::Jt_, std::cout, "Jt");
  // printf("\n");
  return true;
}

template <typename T>
bool LocalTailPosTask<T>::_UpdateTaskJacobian() {
  TK::Jt_ = 0.5 * _robot_sys->_Jc[linkID::HR_abd] +
            0.5 * _robot_sys->_Jc[linkID::HL_abd] -
            0.5 * _robot_sys->_Jc[linkID::HR] -
            0.5 * _robot_sys->_Jc[linkID::HL];

  return true;
}

template <typename T>
bool LocalTailPosTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class LocalTailPosTask<double>;
template class LocalTailPosTask<float>;
