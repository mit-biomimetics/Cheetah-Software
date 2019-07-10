#ifndef JPOS_TASK_Cheetah
#define JPOS_TASK_Cheetah

#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class JPosTask : public Task<T> {
 public:
  JPosTask(const FloatingBaseModel<T>*);
  virtual ~JPosTask();

  DVec<T> _Kp, _Kd;

 protected:
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
