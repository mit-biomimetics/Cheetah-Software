#ifndef KINEMATICS_BODY_POSTURE_TASK
#define KINEMATICS_BODY_POSTURE_TASK

// (Rx, Ry, Rz, X, Y, Z)
#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class BodyPostureTask : public Task<T> {
 public:
  BodyPostureTask(const FloatingBaseModel<T>*);
  virtual ~BodyPostureTask();

  DVec<T> _Kp, _Kd;

 protected:
  // Update op_cmd_
  virtual bool _UpdateCommand(const void* pos_des, const DVec<T>& vel_des,
                              const DVec<T>& acc_des);
  // Update Jt_
  virtual bool _UpdateTaskJacobian();
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot();
  virtual bool _AdditionalUpdate() { return true; }

  const FloatingBaseModel<T>* _robot_sys;
};

#endif
