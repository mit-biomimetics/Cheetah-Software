#ifndef LEG_HEIGHT_TASK
#define LEG_HEIGHT_TASK

// (X, Y, Z)
#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class LegHeightTask : public Task<T> {
 public:
  LegHeightTask(const FloatingBaseModel<T>*, int, int);
  virtual ~LegHeightTask();

  DVec<T> _Kp_kin;
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

  const FloatingBaseModel<T>* _robot_sys;
  int _link_idx;
  int _local_frame_idx;
};

#endif
