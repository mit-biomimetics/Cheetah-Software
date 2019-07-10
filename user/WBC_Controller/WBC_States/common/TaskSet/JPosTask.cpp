#include "JPosTask.hpp"
#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Utilities/Utilities_print.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
JPosTask<T>::JPosTask(const FloatingBaseModel<T>* robot)
    : Task<T>(cheetah::num_act_joint), robot_sys_(robot) {
  TK::Jt_ = DMat<T>::Zero(cheetah::num_act_joint, cheetah::dim_config);
  (TK::Jt_.block(0, 6, cheetah::num_act_joint, cheetah::num_act_joint))
      .setIdentity();
  TK::JtDotQdot_ = DVec<T>::Zero(cheetah::num_act_joint);

  _Kp = DVec<T>::Constant(cheetah::num_act_joint, 50.);
  _Kd = DVec<T>::Constant(cheetah::num_act_joint, 5.);
}

template <typename T>
JPosTask<T>::~JPosTask() {}

template <typename T>
bool JPosTask<T>::_UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                                 const DVec<T>& acc_des) {
  DVec<T>* pos_cmd = (DVec<T>*)pos_des;

  for (size_t i(0); i < cheetah::num_act_joint; ++i) {
    TK::pos_err_[i] = (*pos_cmd)[i] - robot_sys_->_state.q[i];
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * TK::pos_err_[i] +
                     _Kd[i] * (vel_des[i] - robot_sys_->_state.qd[i]) +
                     acc_des[i];
  }
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(op_cmd_, std::cout, "op cmd");
  // pretty_print(*pos_cmd, std::cout, "pos cmd");

  return true;
}

template <typename T>
bool JPosTask<T>::_UpdateTaskJacobian() {
  return true;
}

template <typename T>
bool JPosTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class JPosTask<double>;
template class JPosTask<float>;
