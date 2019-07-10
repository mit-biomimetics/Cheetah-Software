#ifndef STATE_PROVIDER_Cheetah
#define STATE_PROVIDER_Cheetah

#include <cppTypes.h>
#include "Cheetah_DynaCtrl_Definition.h"

template <typename T>
class StateProvider {
 public:
  static StateProvider<T>* getStateProvider();
  ~StateProvider() {}

  int _mode;
  DVec<T> _Q;
  DVec<T> _Qdot;

  size_t _contact_pt[cheetah::num_leg];
  size_t _num_contact;
  Vec3<T> _local_frame_global_pos;

  T _curr_time;
  T _dir_command[2];
  T _ori_command[3];

  int _num_step;

  Vec3<T> _global_body_pos;

  Vec3<T> _global_fr_loc;
  Vec3<T> _global_fl_loc;
  Vec3<T> _global_hr_loc;
  Vec3<T> _global_hl_loc;

 private:
  StateProvider();
};

#endif
