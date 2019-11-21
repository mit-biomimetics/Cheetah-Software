#include "ProblemGenerator.h"
#include "eigenvalues.h"
#include <iostream>

/*!
 * Generate random MPC problem
 */
template<typename T>
QpProblem<T> ProblemGenerator<T>::generateSparseMPC(s64 n_states, s64 n_control, s64 horizon, s64 n_constraint)
{
  std::srand(0);
  QpProblem<T> problem((n_states + n_control) * horizon, (2*n_states + n_control + n_constraint) * horizon);

  DenseMatrix<T> A(n_states, n_states);
  DenseMatrix<T> A1(n_states, n_states);
  DenseMatrix<T> B(n_states, n_control);

  for(s64 attempt = 0; ; attempt++) {
    // step 1. create dynamics matrix A
    A.setIdentity(); // discrete time dynamics, start with I
    A = A + T(0.1) * DenseMatrix<T>::Random(n_states, n_states); // and add some dense random terms everywhere

    // for now, only consider stable systems, so adjust eigenvalues to be magnitude 1 at most.
    A1 = constrainEigenvalueMagnitude<T>(A, T(1));

    // step 2. create control matrix B

    B.setRandom();

    // step 2.5, compute controllability matrix
    DenseMatrix<T> C(n_states, n_states), Ai(n_states, n_states);
    C.setZero();
    Ai.setIdentity();

    for(s64 i = 0; i < n_states - 1; i++) {
      C += Ai * B * B.transpose() * (Ai.transpose());
      Ai = Ai * A1;
    }

    // check if controllable
    T det = C.determinant();
    if(std::abs(det) < 1e-6) {
      printf("Attempt %ld to create dynamics matrix has failed!\n", attempt);
    } else {
      printf("Success on attempt %ld, determinant is %f\n", attempt, det);
      break;
    }
  }

  // step 3. create control and state bounds (time varying)
  Vector<T> uLower(n_control * horizon), uUpper(n_control * horizon), xLower(n_states * horizon), xUpper(n_states * horizon);
  uLower.setRandom();
  uUpper.setRandom();

  for(s64 i = 0; i < n_control * horizon; i++) {
    if(uLower(i) > 0) uLower(i) = -0.1;
    if(uUpper(i) < 0) uUpper(i) = 0.1;
  }

  for(s64 i = 0; i < n_states * horizon; i++) {
    // set the state bounds quite large.
    xLower(i) = -100;
    xUpper(i) = 100;
  }

  // step 3.5: costs
  Vector<T> stateCost(n_states * horizon), controlCost(n_control * horizon);
  stateCost.setRandom();
  controlCost.setRandom();

  for(s64 i = 0; i < n_states * horizon; i++)
    stateCost[i] *= stateCost[i];

  for(s64 i = 0; i < n_control * horizon; i++)
    controlCost[i] *= controlCost[i];

  // set up problem data with NANs to detect if we forget anything
  problem.l.setConstant(NAN);
  problem.u.setConstant(NAN);
  problem.q.setConstant(NAN);
  problem.A.setConstant(0); // i don't touch everything in here
  problem.P.setConstant(0);

  // step 4. constraint matrix 1, dynamics constraints
  // first n_state * horizon rows
  problem.l.block(0, 0, n_states * (horizon - 1), 1).setZero();
  problem.u.block(0, 0, n_states * (horizon - 1), 1).setZero();
  for(s64 i = 0; i < horizon - 1; i++) {
    // i -> i + 1
    // Adyn * xdyn_i + Bdyn * udyn_i - xdyn_i+1 = 0
    s64 startRow = i * n_states;
    s64 xiCol = i * n_states;
    s64 xi1Col = (i + 1) * n_states;
    s64 uiCol = (horizon * n_states) + i * n_control;
    problem.A.block(startRow, xiCol,  n_states, n_states) = -A1;
    problem.A.block(startRow, xi1Col, n_states, n_states).setIdentity();
    problem.A.block(startRow, uiCol,  n_states, n_control) = -B;
  }

  // step 4.5 inital state constraint

  problem.l.block(n_states * (horizon - 1), 0, n_states, 1).setRandom();
  problem.l.block(n_states * (horizon - 1), 0, n_states, 1) *= 50;
  problem.u.block(n_states * (horizon - 1), 0, n_states, 1) = problem.l.block(n_states * (horizon - 1), 0, n_states, 1);
  problem.A.block(n_states * (horizon - 1), 0, n_states, n_states).setIdentity();

  // step 5. constraint matrix 2, state constraints
  // from n_state*horizon to n_state*2*horizon rows
  s64 startRow = horizon * n_states;
  problem.A.block(startRow, 0, n_states * horizon, n_states * horizon).setIdentity();
  problem.l.block(startRow, 0, n_states*horizon, 1) = xLower;
  problem.u.block(startRow, 0, n_states*horizon, 1) = xUpper;

  // step 6. constraint matrix 3, control constraints
  // from 2 * n_state * horizon to 2 * n_state*horizon + n_control * horizon
  startRow += horizon * n_states;
  problem.A.block(startRow, n_states * horizon, n_control * horizon, n_control * horizon).setIdentity();
  problem.l.block(startRow, 0, n_control * horizon, 1) = uLower;
  problem.u.block(startRow, 0, n_control * horizon, 1) = uUpper;

  // step 7: state cost
  problem.P.block(0,0,n_states * horizon, n_states * horizon).diagonal() = stateCost;
  problem.P.block(n_states * horizon, n_states * horizon, n_control * horizon, n_control * horizon).diagonal() = controlCost;

  problem.q.setZero();

//  std::cout << "A matrix: \n" << A << "\n\n";
//  std::cout << "A1 matrix: \n" << A1 << "\n";
//  std::cout << "B matrix: \n" << B << "\n";
//  std::cout << "uMin: " << uLower.transpose() << "\n";
//  std::cout << "uMax: " << uUpper.transpose() << "\n";
//  std::cout << "xMin: " << xLower.transpose() << "\n";
//  std::cout << "xMax: " << xUpper.transpose() << "\n";
//
//  std::cout << "\n\n\n";
//  std::cout << "Problem.A: \n" << problem.A << "\n";
//  std::cout << "Problem.u: \n" << problem.u << "\n";
//  std::cout << "Problem.l: \n" << problem.l << "\n";


  return problem;
}

template class ProblemGenerator<float>;
template class ProblemGenerator<double>;