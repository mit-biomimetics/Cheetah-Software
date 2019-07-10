#ifndef QPSOLVER_QPPROBLEM_H
#define QPSOLVER_QPPROBLEM_H

#include <vector>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Sparse>
#include "types.h"
#include "CholeskyDenseSolver.h"
#include "CholeskySparseSolver.h"

// 0.5 * x'Px + q'x
// l <= Ax <= u

template<typename T>
struct QpProblemSettings {
  s64 maxIterations = 1000000;
  T rho = 6;
  T sigma = 1e-6;
  T alpha = 1.7;

  // numerical hacks
  T infty = 1e10;
  T eqlTol = 1e-10;
  T rhoEqualityScale = 1e3;
  T rhoInfty = 1e-6;

  T terminate = 1e-3;

  void print()
  {
    printf("rho: %f\n"
           "sigma: %f\n"
           "alpha: %f\n", rho, sigma, alpha);
  }
};

enum class ConstraintType {
  INFINITE,
  INEQUALITY,
  EQUALITY
};

template<typename T>
struct ConstraintInfo {
  T rho;
  T invRho;
  ConstraintType type;
};

template<typename T>
class QpProblem
{
public:

    QpProblem(s64 n_, s64 m_, bool print_timings = true)
     : A(m_,n_), P(n_,n_),
    l(m_), u(m_), q(n_), 
    n(n_), m(m_),
    _print(print_timings),
    _kkt(n_ + m_, n_ + m_),
    _cholDenseSolver(print_timings),
    _xzTilde(n_ + m_), _y(m_),
    _x0(n_), _x1(n_), _z0(m_), _z1(m_),   
    _Ar(m_), _Pr(n_), _AtR(n_), 
    _deltaY(m_), _AtDeltaY(n_),
    _deltaX(n_), _PDeltaX(n_), _ADeltaX(m_)
    {
        q.setZero();
        _constraintInfos.resize(m);
      (void)(n_);
      (void)(m_);
    }

    void runFromDense(s64 nIterations = -1, bool sparse = false, bool b_print = true);
    void runFromTriples(s64 nIterations = -1, bool b_print = true);

    Vector<T>& getSolution() { return *_x; }



  // public data
  QpProblemSettings<T> settings;
  DenseMatrix<T> A, P;
  std::vector<SparseTriple<T>> A_triples, P_triples;
  Vector<T> l, u, q;
  s64 n, m;

  ~QpProblem() {
    //printf("done!\n");
  }

private:
  void coldStart();
  void computeConstraintInfos();
  void setupLinearSolverCommon();
  void stepSetup();
  void solveLinearSystem();
  void stepX();
  void stepZ();
  void stepY();
  void setupTriples();
  T calcAndDisplayResidual(bool print);
  T infNorm(const Vector<T>& v);

  bool _print;
  DenseMatrix<T> _kkt;
  std::vector<SparseTriple<T>> _kktTriples;

  CholeskyDenseSolver<T> _cholDenseSolver;
  CholeskySparseSolver<T> _cholSparseSolver;
  Eigen::SparseMatrix<T> Asparse, Psparse;

  Vector<T> _xzTilde, _y;
  Vector<T> _x0, _x1, _z0, _z1;
  Vector<T> *_x, *_z, *_xPrev, *_zPrev;
  Vector<T> _Ar, _Pr, _AtR; // residuals
  Vector<T> _deltaY, _AtDeltaY, _deltaX, _PDeltaX, _ADeltaX; // infeasibilities
  std::vector<ConstraintInfo<T>> _constraintInfos;



  bool _hotStarted = false, _sparse = false;
};


#endif //QPSOLVER_QPPROBLEM_H
