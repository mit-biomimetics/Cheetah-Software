#ifndef B_SPLINE_BASIC
#define B_SPLINE_BASIC

#include <assert.h>
#include <stdio.h>
#include <iostream>

#define SP_IS_EQUAL(x, y) (((x) - (y)) * ((x) - (y)) < 1.e-10)
#define SP_SAFE_DELETE_AR(p) \
  if (p) {                   \
    delete[] p;              \
    (p) = NULL;              \
  }

/*!
 * Basic Bspline  <br/>
 * DIM : Dimension of control points <br/>
 * DEGREE : Derivation is going to be 0 when it is over DEGREE <br/>
 * NUM_MIDDLE : Num middle points (the points except initial and final) <br/>
 * <br/>
 * CONST_LEVEL_INI/FIN : constraint level <br/>
 * 0: position <br/>
 * 1: + velocity <br/>
 * 2: + acceleration <br/>
 * <br/>
 * ******************************  WARNING
 * *********************************************** <br/> NumKnots(DEGREE +
 * NUM_MIDDLE + 2 + CONST_LEVEL_INI + CONST_LEVEL_FIN) >= 2 * (DEGREE + 1) <br/>
 * ******************************  WARNING
 * *********************************************** <br/>
 */

template <typename T, int DIM, int DEGREE, int NUM_MIDDLE, int CONST_LEVEL_INI,
          int CONST_LEVEL_FIN>
class BS_Basic {
 public:
  BS_Basic()
      : NumKnots_(DEGREE + NUM_MIDDLE + 2 + CONST_LEVEL_INI + CONST_LEVEL_FIN +
                  1),
        NumCPs_(NUM_MIDDLE + 2 + CONST_LEVEL_INI + CONST_LEVEL_FIN) {
    for (int i(0); i < NumKnots_; ++i) Knots_[i] = 0.;
    for (int i(0); i < NumCPs_; ++i) {
      for (int j(0); j < DIM; ++j) CPoints_[i][j] = 0.;
    }
    if (NumKnots_ < 2 * (DEGREE + 1)) {
      printf("Invalid setup (num_knots, degree): %d, %d\n", NumKnots_, DEGREE);
    }
  }
  ~BS_Basic() {}

  /*!
   * size of T: DIM * CONST_LEVEL_INI (or CONST_LEVEL_FIN) <br/>
   * ex) if dim:3, const level ini: 3(pos, vel, acc) <br/>
   * ini[0 ~ 2]: pos <br/>
   * ini[3 ~ 5]: vel <br/>
   * ini[6 ~ 8]: acc <br/>
   * @param init : initial point information (vector) <br/>
   * @param fin : finial point information (vector) <br/>
   * @param middle_pt   : middle point information (matrix) <br/>
   * @param fin_time : duration of the spline <br/>
   * @return boolean : success <br/>
   */
  bool SetParam(T* init, T* fin, T** middle_pt, T fin_time) {
    _CalcKnot(fin_time);
    _CalcConstrainedCPoints(init, fin, fin_time);
    _CalcCPoints(middle_pt);

    return true;
  }

  /*!
   * get spline position at the given time <br/>
   * If the input time is before 0, it returns the initial <br/>
   * If the input time is after the final time, it returns the final <br/>
   * @param u : time <br/>
   * @return ret : position at the given time.  <br/>
   */
  bool getCurvePoint(T u, T* ret) {
    int _span;

    if (u < Knots_[0])
      u = Knots_[0];
    else if (u > Knots_[NumKnots_ - 1]) {
      u = Knots_[NumKnots_ - 1];
    }

    if (!_findSpan(_span, u)) return false;

    T _N[DEGREE + 1];
    _BasisFuns(_N, _span, u);

    T _C[DIM];

    for (int j(0); j < DIM; ++j) {
      _C[j] = 0.0;
      for (int i(0); i <= DEGREE; ++i) {
        _C[j] += _N[i] * CPoints_[_span - DEGREE + i][j];
      }
    }

    for (int i(0); i < DIM; ++i) {
      ret[i] = _C[i];
    }
    return true;
  }

  /*!
   * get spline derivative information at the given time <br/>
   * If the input time is before 0, it returns the initial <br/>
   * If the input time is after the final time, it returns the final <br/>
   * @param u : time <br/>
   * @param d : drivative level (e.g. 1: velocity, 2: acceleration) <br/>
   * @return ret : derivative information at the given time.  <br/>
   */
  bool getCurveDerPoint(T u, int d, T* ret) {
    if (d > DEGREE) return 0.0;

    if (u < Knots_[0])
      u = Knots_[0];
    else if (u > Knots_[NumKnots_ - 1]) {
      u = Knots_[NumKnots_ - 1];
    }

    T** _CK = new T*[d + 1];
    for (int i(0); i < d + 1; ++i) {
      _CK[i] = new T[DIM];
    }
    if (_CurveDerivsAlg1V(_CK, u, d)) {
      for (int m(0); m < DIM; ++m) ret[m] = _CK[d][m];

      for (int p(0); p < d + 1; ++p) delete[] _CK[p];
      SP_SAFE_DELETE_AR(_CK);
      return true;
    } else {
      for (int p(0); p < d + 1; ++p) delete[] _CK[p];
      SP_SAFE_DELETE_AR(_CK);
    }
    return false;
  }

  // protected:
 private:
  inline void _CalcKnot(T Tf) {
    int _i(0);
    int _j(0);
    int _NumMidKnot(NumKnots_ - 2 * DEGREE - 2);
    T _TimeStep = Tf / (_NumMidKnot + 1);

    // augment knot sequence for the initial part, # of order ( degree + 1 )
    for (_j = 0; _j < DEGREE + 1; ++_j) Knots_[_i++] = 0.0;

    // uniform knot sequence for the middle part,
    // #: NumKnot - degree - degree = NumKnot - order - order + 2
    for (_j = 0; _j < _NumMidKnot; ++_j) {
      Knots_[_i] = Knots_[_i - 1] + _TimeStep;
      ++_i;
    }
    // augment knot sequence for the final part, # of order ( degree + 1 )
    for (_j = 0; _j < DEGREE + 1; ++_j) Knots_[_i++] = Tf;

    // for(int i(0); i< NumKnots_; ++i)
    // std::cout<<Knots_[i]<<std::endl;
  }

  bool _CurveDerivsAlg1V(T** CK, T u, int d) {
    assert(d <= DEGREE);

    int _k(0);
    int _j(0);

    T** _nders = new T*[d + 1];

    for (_k = 0; _k < d + 1; ++_k) _nders[_k] = new T[DEGREE + 1];

    int _span;
    if (!_findSpan(_span, u)) return false;

    _BasisFunsDers(_nders, _span, u, d);

    for (_k = 0; _k <= d; ++_k) {
      // Clean Up Column
      for (int m(0); m < DIM; ++m) CK[_k][m] = 0.;

      for (_j = 0; _j <= DEGREE; ++_j) {
        for (int m(0); m < DIM; ++m) {
          CK[_k][m] += _nders[_k][_j] * CPoints_[_span - DEGREE + _j][m];
        }
      }
    }
    for (_k = 0; _k < d + 1; ++_k) delete[] _nders[_k];
    delete[] _nders;

    return true;
  }
  bool _BasisFunsDers(T** ders, T u, int n) {
    // int _span = FindSpan(u);
    int _span;
    if (!_findSpan(_span, u)) return false;

    _BasisFunsDers(ders, _span, u, n);
    return true;
  }

  bool _BasisFunsDers(T** ders, int span, T u, int n) {
    int _j, _r, _k;
    int _s1, _s2;
    int _j1, _j2;
    int _rk;
    int _pk;

    T _saved = 0.0;
    T _left = 0.0;
    T _right = 0.0;
    T _temp = 0.0;
    T _d = 0.0;

    // to store the basis functions and knot differences
    T** _ndu = new T*[DEGREE + 1];
    for (_j = 0; _j <= DEGREE; ++_j) _ndu[_j] = new T[DEGREE + 1];
    // to store (in an alternating fashion) the two most recently computed
    // rows a(k,j) and a(k-1,j)
    T** _a = new T*[2];
    for (_j = 0; _j < 2; ++_j) _a[_j] = new T[DEGREE + 1];

    _ndu[0][0] = 1.0;
    for (_j = 1; _j <= DEGREE; ++_j) {
      _saved = 0.0;
      for (_r = 0; _r < _j; ++_r) {
        _left = _Left(span, _j - _r, u);
        _right = _Right(span, _r + 1, u);

        // Lower triangle
        _ndu[_j][_r] = _right + _left;
        _temp = _ndu[_r][_j - 1] / _ndu[_j][_r];

        // Upper triangle
        _ndu[_r][_j] = _saved + _right * _temp;
        _saved = _left * _temp;
      }
      _ndu[_j][_j] = _saved;
    }

    // Load the basis functions
    for (_j = 0; _j <= DEGREE; ++_j) ders[0][_j] = _ndu[_j][DEGREE];

    // This section computes the derivatives (Eq. [2.9])
    for (_r = 0; _r <= DEGREE; ++_r) {
      _s1 = 0;
      _s2 = 1;
      _a[0][0] = 1.0;

      // Loop to compute k-th derivative
      for (_k = 1; _k <= n; ++_k) {
        _d = 0.0;
        _rk = _r - _k;
        _pk = DEGREE - _k;

        if (_r >= _k) {
          _a[_s2][0] = _a[_s1][0] / _ndu[_pk + 1][_rk];
          _d = _a[_s2][0] * _ndu[_rk][_pk];
        }

        if (_rk >= -1)
          _j1 = 1;
        else
          _j1 = -_rk;

        if (_r - 1 <= _pk)
          _j2 = _k - 1;
        else
          _j2 = DEGREE - _r;

        for (_j = _j1; _j <= _j2; ++_j) {
          _a[_s2][_j] =
              (_a[_s1][_j] - _a[_s1][_j - 1]) / _ndu[_pk + 1][_rk + _j];
          _d += _a[_s2][_j] * _ndu[_rk + _j][_pk];
        }

        if (_r <= _pk) {
          _a[_s2][_k] = -_a[_s1][_k - 1] / _ndu[_pk + 1][_r];
          _d += _a[_s2][_k] * _ndu[_r][_pk];
        }
        ders[_k][_r] = _d;

        // Switch rows
        _j = _s1;
        _s1 = _s2;
        _s2 = _j;
      }
    }

    // Multiply through by the correct factors
    // (Eq. [2.9])
    _r = DEGREE;
    for (_k = 1; _k <= n; ++_k) {
      for (_j = 0; _j <= DEGREE; ++_j) ders[_k][_j] *= _r;

      _r *= (DEGREE - _k);
    }

    // Deallocate
    for (_j = 0; _j <= DEGREE; ++_j) delete[] _ndu[_j];
    delete[] _ndu;

    for (_j = 0; _j < 2; ++_j) delete[] _a[_j];
    delete[] _a;

    return true;
  }

  void _BasisFuns(T* N, T u) {
    // Original
    /*int _span = FindSpan(u);
      BasisFuns(N, _span, u);*/

    int _span;

    if (_findSpan(_span, u)) {
      _BasisFuns(N, _span, u);
    }
  }

  void _BasisFuns(T* N, int span, T u) {
    int _j, _r;
    T _left = 0.0;
    T _right = 0.0;
    T _saved = 0.0;
    T _temp = 0.0;

    N[0] = 1.0;
    for (_j = 1; _j <= DEGREE; ++_j) {
      _saved = 0.0;
      for (_r = 0; _r < _j; ++_r) {
        _left = _Left(span, _j - _r, u);
        _right = _Right(span, _r + 1, u);

        if ((_right + _left) != 0) {
          _temp = N[_r] / (_right + _left);
        }

        N[_r] = _saved + _right * _temp;
        _saved = _left * _temp;
      }
      N[_j] = _saved;
    }
  }
  inline T _Left(int i, int j, T u) { return u - Knots_[i + 1 - j]; }

  inline T _Right(int i, int j, T u) { return Knots_[i + j] - u; }

  bool _findSpan(int& ret, T u) {
    if (u < Knots_[0] || Knots_[NumKnots_ - 1] < u) return false;

    if (SP_IS_EQUAL(u, Knots_[NumKnots_ - 1])) {
      for (int i(NumKnots_ - 2); i > -1; --i) {
        if (Knots_[i] < u && u <= Knots_[i + 1]) {
          ret = i;
          return true;
        }
      }
      return false;
    }
    // Binary search
    int _low = 0;
    int _high = NumKnots_ - 1;
    int _mid = (_low + _high) >> 1;

    while (u < Knots_[_mid] || u >= Knots_[_mid + 1]) {
      if (u < Knots_[_mid])
        _high = _mid;
      else
        _low = _mid;
      _mid = (_low + _high) >> 1;
    }
    ret = _mid;
    return true;
  }

  void _CalcConstrainedCPoints(T* init, T* fin, T Tf) {
    // Position
    for (int m(0); m < DIM; ++m) {
      CPoints_[0][m] = init[m];
      CPoints_[NumCPs_ - 1][m] = fin[m];
    }
    // Initial Constraints
    T** d_mat = new T*[CONST_LEVEL_INI + 1];

    for (int i(0); i < CONST_LEVEL_INI + 1; ++i)
      d_mat[i] = new T[CONST_LEVEL_INI + 2];

    _BasisFunsDers(d_mat, 0., CONST_LEVEL_INI);

    T ini_const[DIM];
    // Vel, Acc, ...
    for (int j(1); j < CONST_LEVEL_INI + 1; ++j) {
      for (int k(0); k < DIM; ++k) {
        ini_const[k] = init[j * DIM + k];

        for (int h(j); h > 0; --h) {
          ini_const[k] -= d_mat[j][h - 1] * CPoints_[h - 1][k];
        }
        CPoints_[j][k] = ini_const[k] / d_mat[j][j];
      }
    }

    for (int p(0); p < CONST_LEVEL_INI + 1; ++p) delete[] d_mat[p];
    SP_SAFE_DELETE_AR(d_mat);

    // Final Constraints
    T** c_mat = new T*[CONST_LEVEL_FIN + 1];

    for (int i(0); i < CONST_LEVEL_FIN + 1; ++i)
      c_mat[i] = new T[CONST_LEVEL_FIN + 2];

    _BasisFunsDers(c_mat, Tf, CONST_LEVEL_FIN);

    // Vel, Acc, ...
    int idx(1);
    for (int j(NumCPs_ - 2); j > NumCPs_ - 2 - CONST_LEVEL_FIN; --j) {
      for (int k(0); k < DIM; ++k) {
        ini_const[k] = fin[idx * DIM + k];

        for (int h(idx); h > 0; --h) {
          ini_const[k] -=
              c_mat[idx][CONST_LEVEL_FIN + 2 - h] * CPoints_[NumCPs_ - h][k];
        }
        CPoints_[j][k] = ini_const[k] / c_mat[idx][CONST_LEVEL_FIN + 1 - idx];
      }
      ++idx;
    }
    for (int p(0); p < CONST_LEVEL_FIN + 1; ++p) delete[] c_mat[p];
    SP_SAFE_DELETE_AR(c_mat);
  }

  void _PrintCP(int i) {
    printf("%i th CP check:\n", i);
    for (int m(0); m < DIM; ++m) {
      for (int j(0); j < NumCPs_; ++j) {
        printf("%f \t", CPoints_[j][m]);
      }
      printf("\n");
    }
    printf("\n");
  }

  void _CalcCPoints(T** middle_pt) {
    for (int i(0); i < NUM_MIDDLE; ++i) {
      for (int m(0); m < DIM; ++m) {
        CPoints_[CONST_LEVEL_INI + 1 + i][m] = middle_pt[i][m];
      }
    }
  }

  T fin_time_;
  int NumKnots_;
  int NumCPs_;

  T Knots_[DEGREE + NUM_MIDDLE + 2 + CONST_LEVEL_INI + CONST_LEVEL_FIN + 1];
  T CPoints_[NUM_MIDDLE + 2 + CONST_LEVEL_INI + CONST_LEVEL_FIN][DIM];
};

#endif
