#include "ImpulseCurve.hpp"

template <typename T>
ImpulseCurve<T>::ImpulseCurve() {}
template <typename T>
ImpulseCurve<T>::~ImpulseCurve() {}

template <typename T>
void ImpulseCurve<T>::setCurve(const T& apex, const T& time) {
  _end_time = time;
  T** ctrl_pt = new T*[num_ctrl_pt];
  for (int i(0); i < num_ctrl_pt; ++i) {
    ctrl_pt[i] = new T[dim];
  }

  ctrl_pt[0][0] = 0.;
  ctrl_pt[1][0] = 0.8 * apex;
  ctrl_pt[2][0] = 1.0 * apex;
  ctrl_pt[3][0] = 1.0 * apex;
  _increasing_part.SetParam(ctrl_pt, 0.5);

  ctrl_pt[0][0] = 1.0 * apex;
  ctrl_pt[1][0] = 1.0 * apex;
  ctrl_pt[2][0] = 0.8 * apex;
  ctrl_pt[3][0] = 0.;
  _decreasing_part.SetParam(ctrl_pt, 0.5);

  for (int i(0); i < num_ctrl_pt; ++i) {
    delete[] ctrl_pt[i];
  }
  delete[] ctrl_pt;
}

template <typename T>
T ImpulseCurve<T>::getValue(const T& t) {
  T ret[dim];
  T u = t / _end_time;

  if (t < _end_time / 2.) {
    _increasing_part.getCurvePoint(u, ret);
  } else {
    _decreasing_part.getCurvePoint(u - 0.5, ret);
  }
  return ret[0];
}

template class ImpulseCurve<double>;
template class ImpulseCurve<float>;
