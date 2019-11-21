#ifndef IMPULSE_CURVE
#define IMPULSE_CURVE

#include <Utilities/BezierCurve.h>

template <typename T>
class ImpulseCurve {
 public:
  ImpulseCurve();
  ~ImpulseCurve();

  void setCurve(const T& apex, const T& time);
  T getValue(const T& t);

 private:
  constexpr static int num_ctrl_pt = 4;
  constexpr static int dim = 1;
  T _end_time;
  T _apex_value;

  BezierCurve<T, dim, num_ctrl_pt> _increasing_part;
  BezierCurve<T, dim, num_ctrl_pt> _decreasing_part;
};
#endif
