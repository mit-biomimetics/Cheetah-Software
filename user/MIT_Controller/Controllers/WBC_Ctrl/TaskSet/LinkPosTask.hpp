#ifndef LINK_POS_TASK
#define LINK_POS_TASK

// (X, Y, Z)
#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class LinkPosTask : public Task<T> {
 public:
  LinkPosTask(const FloatingBaseModel<T>*, int link_idx,
              bool virtual_depend = true);
  virtual ~LinkPosTask();

  DVec<T> _Kp, _Kd, _Kp_kin;

 protected:
  // Update op_cmd_
  virtual bool _UpdateCommand(const void* pos_des, const DVec<T>& vel_des,
                              const DVec<T>& acc_des);
  // Update Jt_
  virtual bool _UpdateTaskJacobian();
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot();
  virtual bool _AdditionalUpdate() { return true; }

  const FloatingBaseModel<T>* robot_sys_;
  int link_idx_;
  bool virtual_depend_;
};

#endif
