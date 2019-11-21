#ifndef BEZIER_CURVE
#define BEZIER_CURVE

#include <math.h>

// N = NUM_CTR_PT -1
template <typename T, int DIM, int NUM_CTR_PT>
class BezierCurve {
 public:
  T _CtrlPt[NUM_CTR_PT][DIM];
  T _coeff[NUM_CTR_PT];
  T _end_time;
  BezierCurve() {
    for (int j(0); j < NUM_CTR_PT; ++j) {
      for (int i(0); i < DIM; ++i) {
        _CtrlPt[j][i] = 0.;
      }
      _coeff[j] = 0.;
    }
  }

  // ctrl pt 0: initial
  // ctrl pt n-1(end): final
  bool SetParam(T** ctrl_pt, T fin_time) {
    _end_time = fin_time;
    T n_fact = _factorial(NUM_CTR_PT - 1);
    for (int j(0); j < NUM_CTR_PT; ++j) {
      for (int i(0); i < DIM; ++i) {
        _CtrlPt[j][i] = ctrl_pt[j][i];
      }
      _coeff[j] = n_fact / (_factorial(j) * _factorial((NUM_CTR_PT - 1) - j));
    }
    return true;
  }

  bool getCurvePoint(T u, T* ret) {
    if (u > _end_time) {
      for (int i(0); i < DIM; ++i) {
        ret[i] = _CtrlPt[NUM_CTR_PT - 1][i];
      }
      return true;
    } else if (u < 0.) {
      for (int i(0); i < DIM; ++i) {
        ret[i] = _CtrlPt[0][i];
      }
      return true;
    } else {
      u /= _end_time;
      for (int i(0); i < DIM; ++i) {
        ret[i] = 0.;
        for (int j(0); j < NUM_CTR_PT; ++j) {
          ret[i] += _coeff[j] * pow(u, j) * pow(1 - u, (NUM_CTR_PT - 1) - j) *
                    _CtrlPt[j][i];
        }
      }
    }
    return true;
  }

  bool getCurveVelocity(T u, T* ret) {
    if (u > _end_time) {
      for (int i(0); i < DIM; ++i) {
        ret[i] = 0.;
      }
      return true;
    } else if (u < 0.) {
      for (int i(0); i < DIM; ++i) {
        ret[i] = 0.;
      }
      return true;
    } else {
      u /= _end_time;
      for (int i(0); i < DIM; ++i) {
        ret[i] = 0.;
        ret[i] += _coeff[0] *
                  (-(NUM_CTR_PT - 1) * pow(1 - u, (NUM_CTR_PT - 1) - 1)) *
                  _CtrlPt[0][i];
        for (int j(1); j < NUM_CTR_PT - 1; ++j) {
          ret[i] += _coeff[j] *
                    (j * pow(u, j - 1) * pow(1 - u, (NUM_CTR_PT - 1) - j) -
                     (NUM_CTR_PT - 1 - j) * pow(u, j) *
                         pow(1 - u, (NUM_CTR_PT - 1) - j - 1)) *
                    _CtrlPt[j][i];
        }
        ret[i] += _coeff[NUM_CTR_PT - 1] * (NUM_CTR_PT - 1) *
                  pow(u, NUM_CTR_PT - 2) * _CtrlPt[NUM_CTR_PT - 1][i];

        ret[i] /= _end_time;
      }
    }
    return true;
  }

 private:
  int _factorial(const int& n) {
    int ret = 1;
    for (int i(1); i < n + 1; ++i) ret *= (i);

    return ret;
  }
};
#endif
