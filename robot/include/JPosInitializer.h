#ifndef JPOS_INITIALIZER
#define JPOS_INITIALIZER

#include <Controllers/LegController.h>
#include <Dynamics/Quadruped.h>
#include <Utilities/BSplineBasic.h>

template <typename T>
class JPosInitializer {
 public:
  JPosInitializer(T end_time);
  ~JPosInitializer();

  bool IsInitialized(LegController<T>*);

 private:
  void _UpdateParam();
  void _UpdateInitial(const LegController<T>* ctrl);
  bool _b_first_visit;
  T _end_time;
  T _curr_time;
  T _dt;

  std::vector<T> _ini_jpos;
  std::vector<T> _target_jpos;
  std::vector<T> _mid_jpos;

  BS_Basic<T, cheetah::num_act_joint, 3, 1, 2, 2> _jpos_trj;
};
#endif
