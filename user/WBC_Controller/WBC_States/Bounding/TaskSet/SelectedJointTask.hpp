#ifndef SELECTED_JPOS_TASK_Cheetah
#define SELECTED_JPOS_TASK_Cheetah

#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class SelectedJointTask : public Task<T> {
 public:
  SelectedJointTask(const FloatingBaseModel<T>*, const std::vector<int>&);
  virtual ~SelectedJointTask();

  DVec<T> _Kp, _Kd;

 protected:
  std::vector<int> _jidx_list;
  // Update op_cmd_
  virtual bool _UpdateCommand(void* pos_des, const DVec<T>& vel_des,
                              const DVec<T>& acc_des);
  // Update Jt_
  virtual bool _UpdateTaskJacobian();
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot();
  virtual bool _AdditionalUpdate() { return true; }
  const FloatingBaseModel<T>* robot_sys_;
};

#endif
