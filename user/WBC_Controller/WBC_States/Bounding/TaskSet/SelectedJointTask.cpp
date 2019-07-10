#include "SelectedJointTask.hpp"
#include <Configuration.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Utilities/Utilities_print.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
SelectedJointTask<T>::SelectedJointTask(const FloatingBaseModel<T>* robot,
                                        const std::vector<int>& jidx_list)
    : Task<T>(jidx_list.size()), _jidx_list(jidx_list), robot_sys_(robot) {
  TK::Jt_ = DMat<T>::Zero(TK::dim_task_, cheetah::dim_config);
  for (size_t i(0); i < TK::dim_task_; ++i) TK::Jt_(i, _jidx_list[i] + 6) = 1.;
  TK::JtDotQdot_ = DVec<T>::Zero(TK::dim_task_);

  _Kp = DVec<T>::Constant(TK::dim_task_, 70.);
  _Kd = DVec<T>::Constant(TK::dim_task_, 5.);

  // printf("dim: %lu\n", TK::dim_task_);
}

template <typename T>
SelectedJointTask<T>::~SelectedJointTask() {}

template <typename T>
bool SelectedJointTask<T>::_UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                                          const DVec<T>& acc_des) {
  DVec<T>* pos_cmd = (DVec<T>*)pos_des;

  for (size_t i(0); i < TK::dim_task_; ++i) {
    TK::pos_err_[i] = (*pos_cmd)[i] - robot_sys_->_state.q[_jidx_list[i]];
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] =
        _Kp[i] * TK::pos_err_[i] +
        _Kd[i] * (TK::vel_des_[i] - robot_sys_->_state.qd[_jidx_list[i]]) +
        TK::acc_des_[i];
  }
  // pretty_print(*pos_cmd, std::cout, "pos cmd");
  // pretty_print(vel_des, std::cout, "vel_des");
  // pretty_print(acc_des, std::cout, "acc_des");
  // pretty_print(TK::op_cmd_, std::cout, "op cmd");
  // pretty_print(TK::pos_err_, std::cout, "pos err");
  // pretty_print(robot_sys_->_state.q, std::cout, "jpos");
  // pretty_print(robot_sys_->_state.qd, std::cout, "jvel");
  // printf("\n");
  return true;
}

template <typename T>
bool SelectedJointTask<T>::_UpdateTaskJacobian() {
  return true;
}

template <typename T>
bool SelectedJointTask<T>::_UpdateTaskJDotQdot() {
  return true;
}

template class SelectedJointTask<double>;
template class SelectedJointTask<float>;
