#include "JumpNLP.hpp"
#include <iostream>

#ifdef IPOPT_OPTION

#include <cassert>

using namespace Ipopt;

template<typename T>
JumpNLP<T>::JumpNLP(){}


template <typename T>
bool JumpNLP<T>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
    Index& nnz_h_lag, IndexStyleEnum& index_style){
  // Number of total decision variables for the problem
  n = NUM_DECISION_VARS * NUM_PREDICTIONS;

  // Number of total constraints for the problem
  m = NUM_CONSTRAINTS;

  // Number of nonzeros in the constraint Jacobian
  nnz_jac_g = NNZ_J_G * (NUM_PREDICTIONS - 1) + NNZ_J_G_F; // Last step only dependent on current so dynamics goes away

  // Number of nonzeros in the Hessian matrix (number calculated from MATLAB)
  nnz_h_lag = NNZ_H * (NUM_PREDICTIONS - 1) + NNZ_H_F; // Currently only for cost hessian, need constraints
  // Will definitely change possibly... It didnt change! Maybe...

  // Index style is 0-based
  index_style = C_STYLE;

  // Set the new variables from the user input
  SetParameters();

  return true;
}

/** Method to return the bounds for my problem */
template <typename T>
bool JumpNLP<T>::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u){
  assert(n == n_opt);
  assert(m == n_constr);

  Bounds(x_u, x_l, g_u, g_l);
  return true;
}

/** Method to return the starting point for the algorithm */
template <typename T>
bool JumpNLP<T>::get_starting_point(Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,
    Number* lambda){
  (void)n;
  (void)z_L;
  (void)z_U;
  (void)m;
  (void)lambda;

  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  InitializeMPC(x);

  return true;
}

/** Method to return the objective value */
template <typename T>
bool JumpNLP<T>::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
  //bool JumpNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
  (void)n;
  (void)new_x;

  obj_value = PredictedCost(x);

  return true;
}


/** Method to return the gradient of the objective */
template <typename T>
bool JumpNLP<T>::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
  (void)n;
  (void)new_x;
  PredictedCostGradient(x, grad_f);
  return true;
}

/** Method to return the constraint residuals */
template <typename T>
bool JumpNLP<T>::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
  (void)n;
  (void)new_x;
  (void)m;

  Constraints(x, g);

  return true;
}

/** Method to return:
 *   1) The structure of the jacobian (if "values" is NULL)
 *   2) The values of the jacobian (if "values" is not NULL)
 */
template <typename T>
bool JumpNLP<T>::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol,
    Number* values){
  (void)n;
  (void)x;
  (void)new_x;
  (void)m;
  (void)nele_jac;

  if (values == NULL) {
    //cout << "[PRMPC] Constraint Jacobian Setup Start" << endl;

    for (int i = 0; i < NUM_PREDICTIONS; i++) {
      i_jc = i * NNZ_J_G;
      if (i < (NUM_PREDICTIONS - 1)) {
        Jump2DConstraintJacobianSP(
            i, NUM_DECISION_VARS, NUM_CONSTRAINTS, iRow + i_jc, jCol + i_jc);
      } else {
        Jump2DConstraintJacobianFinalSP(
            i, NUM_DECISION_VARS, NUM_CONSTRAINTS, iRow + i_jc, jCol + i_jc);
      }
    }
    //cout << "[PRMPC] Constraint Jacobian Setup End" << endl;

  } else {
    //cout << "[PRMPC] Constraint Jacobian Calculation Start" << endl;

    // Calculate the constraint jacobian
    ConstraintJacobian(x, values);
    //cout << "[PRMPC] Constraint Jacobian Calculation End" << endl;

  }


  return true;
}

/** Method to return:
 *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
 *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
 */
template <typename T>
bool JumpNLP<T>::eval_h(Index n, const Number* x, bool new_x,
    Number obj_factor, Index m, const Number* lambda,
    bool new_lambda, Index nele_hess, Index* iRow,
    Index* jCol, Number* values){
  (void)n;
  (void)x;
  (void)new_x;
  (void)m;
  (void)new_lambda;
  (void)nele_hess;

  if (values == NULL) {
    //cout << "[PRMPC] Hessian Setup Start" << endl;

    for (int i = 0; i < NUM_PREDICTIONS; i++) {
      i_h = i * NNZ_H;

      if (i < (NUM_PREDICTIONS - 1)) {
        Jump2DLagrangianHessianSP(i, NUM_DECISION_VARS, iRow + i_h, jCol + i_h);
      } else {
        Jump2DLagrangianHessianFinalSP(i, NUM_DECISION_VARS, iRow + i_h, jCol + i_h);
      }
    }
    //cout << "[PRMPC] Hessian Setup End" << endl;

  } else {
    //cout << "[PRMPC] Hessian Calculation Start" << endl;

    // Calculate the Hessian
    LagrangianHessian(x, values, obj_factor, lambda);
    //cout << "[PRMPC] Hessian Calculation End" << endl;

  }

  return true;
}

/** @name Solution Methods */
/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
template <typename T>
void JumpNLP<T>::finalize_solution(SolverReturn status,
    //void JumpNLP::finalize_solution(SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda,
    Number obj_value, const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq){

    (void)status;
    (void)n;
    (void)x;
    (void)z_L;
    (void)z_U;
    (void)m;
    (void)g;
    (void)lambda;
    (void)obj_value;
    (void)ip_data;
    (void)ip_cq;

    printf("solution: %f, %f\n", x[0], x[1]);
    printf("solver return: %d\n", status);
}

/**
 * Create bounds for the decision variables and the constraints.
 * Note that the first timestep states are constrained to the robot's
 * state at the beginning of the optimization.
 */
template <typename T>
void JumpNLP<T>::Bounds(double* X_ub, double* X_lb, double* constraints_ub, double* constraints_lb) {

  // Create the bounds for the prediction horizon
  for (int k = 0; k < NUM_PREDICTIONS; k++) {
    i_xd = k * NUM_STATES; // eventually get an input trajectory
    i_X = k * NUM_DECISION_VARS;
    i_c = k * NUM_CONSTRAINTS;
    i_f = k * NUM_FEET;

    // First step bounds have to be current state
    if (k == 0) {
      Jump2DBounds(
          current_state, current_state, inputs_max, inputs_min, contact_state_pred + i_f,
          X_lb, X_ub, constraints_ub, constraints_lb);

    } else {
      Jump2DBounds(
          states_max, states_min, inputs_max, inputs_min, contact_state_pred + i_f,
          X_lb + i_X, X_ub + i_X, constraints_ub + i_c, constraints_lb + i_c);
    }
  }
}


/**
 * Sets the seeded state trajectory, x_seed, and the seeded input reference
 * policy, u_seed, by calculating the reference policy, (r_ref, f_ref), for
 * all of the feet throughout the prediction horizon trajectory.
 *
 * TODO:
 *  - Add turning feed forward
 */
template <typename T>
void JumpNLP<T>::InitializeMPC(double* X_0) {

  // Set up a temporary vector to hold the current output
  double x_0[NUM_DECISION_VARS] = {0};



  // Copy the current state into the temporary variable
  for (int i = 0; i < NUM_STATES; i++) {
    x_0[i] = current_state[i];
  }


  //cout << contact_state_pred[0] << " " << " " << touchdown_pred[0] << " "  << T_stance[0] << " " << x_d[6] << endl;
  // Create the initial guess and policy seeding
  for (int k = 0; k < NUM_PREDICTIONS; k++) {

    // Start indices
    i_x = k * NUM_DECISION_VARS;
    i_x1 = (k + 1) * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_xd = k * NUM_STATES; // eventually get an input trajectory
    i_uref = k * NUM_INPUTS;
    i_f = k * NUM_FEET;

    /*for (int i = 0; i < NUM_FEET; i++) {

      cout << _gravity << " " << T_stance[0] <<" " << touchdown_pred[i_f+i]<<" " <<p_foot_0[i*3+0] << " " << p_foot_0[i*3+1] << " " << p_foot_0[i*3+2] << endl;;
      }*/
    // Set the current decision variables for the states
    for (int i = 0; i < NUM_STATES; i++) {
      X_0[i_x + i] = x_0[i];
    }

    // Mark the location at the beginning of touchdown
    for (int foot = 0; foot < NUM_FEET; foot++) {
      if (touchdown_pred[i_f + foot] == 1 || k == 0) {
        for (int i = 0; i < 2; i++) {
          x_touchdown[foot * 2 + i] = x_0[i]; // CHECK THIS... Prob in testing
        }
      }
    }

    // Iterative initialization function
    //Jump2DInitialize();
    /*for (int i = 0; i < NUM_FEET; i++) {
      cout << p_foot[i*3+0] << " " << p_foot_0[i*3+1] << " " << p_foot_0[i*3+2];
      cout << endl;
      }*/
    /*
       for (int i = 0; i < NUM_FEET; i++) {

       cout << _gravity << " " << T_stance[0] <<" " << touchdown_pred[i_f+i]<<" " <<p_foot_0[i*3+0] << " " << p_foot_0[i*3+1] << " " << p_foot_0[i*3+2] << endl;;
       }
       cout<< endl;*/

    // Set the current decision variables for the inputs
    for (int i = 0; i < NUM_INPUTS; i++) {
      X_0[i_u + i] = x_0[NUM_STATES + i];
    }
  }
  /*
     cout << "initializing" << endl;
     for (int ind = 0; ind < NUM_INPUTS; ind++){
     for (int ind2 = 0; ind2 < NUM_PREDICTIONS; ind2++){
     cout << std::right << std::setw(13) << u_ref[ind + ind2*NUM_INPUTS];
     }
     cout << "" << endl;
     }
     cout << "\nX0" << endl;
     for (int ind = 0; ind < NUM_DECISION_VARS; ind++){
     for (int ind2 = 0; ind2 < NUM_PREDICTIONS; ind2++){
     cout << std::right << std::setw(13) << X_0[ind + ind2*NUM_DECISION_VARS];
     }
     cout << "" << endl;
     }
     */
  /*
     cout << "\nTD" << endl;
     for (int ind = 0; ind < NUM_FEET; ind++){
     for (int ind2 = 0; ind2 < NUM_PREDICTIONS; ind2++){
     cout << std::right << std::setw(13) << touchdown_pred[ind + ind2*NUM_FEET];
     }
     cout << "" << endl;
     }
     cout << "\n" << endl;*/
}

/**
 * Predicted cost over the given gait-generalized prediction horizon.
 * Calls the MPCCost function that is generated first symbolically in
 * MATLAB to produce an optimized MATLAB function, then a function is
 * created using MATLAB's codegen. This method is essentially a wrapper
 * that takes the problem data and uses this iterative cost function with
 * the known gait logic.
 *
 *
 *  J(X) = x_error'*Q*x_error + u_error'*R*u_error
 *
 *      x_error = x_d - x
 *      u_error = u_ref - u
 */
template<typename T>
double JumpNLP<T>::PredictedCost(const double * X) {

  // Initialize the predicted cost
  double cost = 0;

  // Incremental cost over the predicted gait horizon
  for (int k = 0; k < NUM_PREDICTIONS; k++) {

    // Start indices
    i_x = k * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_xd = k * NUM_STATES; // eventually get an input trajectory
    i_uref = k * NUM_INPUTS;
    i_f = k * NUM_FEET;
    i_Q = k * NUM_STATES;
    i_R = k * NUM_INPUTS;


    //cout << x_d[2] << endl;

    // Calculate the incremental cost
    cost = cost + Jump2DCost(X + i_x, X + i_u, x_d + i_xd, u_ref + i_uref, Q + i_Q, R + i_R);
  }
  /*
     for (int ind = 0; ind < NUM_DECISION_VARS; ind++){
     for (int ind2 = 0; ind2 < NUM_PREDICTIONS; ind2++){
     cout << std::right << std::setw(13) << X[ind + ind2*NUM_DECISION_VARS];
     }
     cout << "" << endl;
     }*/

  //std::usleep(10000000);

  // Return the cost over the prediction horizon
  return cost;
}


/**
 * Predicted cost gradient over the given gait-generalized prediction horizon.
 */
template<typename T>
void JumpNLP<T>::PredictedCostGradient(const double * X, double * gradient) {

  // Incremental gradient over the predicted gait horizon
  for (int k = 0; k < NUM_PREDICTIONS; k++) {
    // Start indices
    i_x = k * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_xd = k * NUM_STATES;
    i_uref = k * NUM_INPUTS;
    i_g = k * NUM_DECISION_VARS;
    i_f = k * NUM_FEET;
    i_Q = k * NUM_STATES;
    i_R = k * NUM_INPUTS;

    // Calculate the iterative gradient of the cost w.r.t. the decision variables
    Jump2DCostGradient(X + i_x, X + i_u, x_d + i_xd, u_ref + i_uref,
        Q + i_Q, R + i_R,
        gradient + i_g);

  }
}

template<typename T>
void JumpNLP<T>::SetParameters() {
}

template<typename T>
void JumpNLP<T>::SetFootLocations(double * p_foot_0_in) {
  //pthread_mutex_lock(&solve_mutex);
  for (int foot = 0; foot < NUM_FEET; foot++) {
    for (int i = 0; i < 2; i++) {
      if (i == 0 || i == 1) {
        p_foot_0[foot * 2 + i] = p_foot_0_in[foot * 2 + i]; // - current_state[i];
      } else {
        p_foot_0[foot * 2 + i] = p_foot_0_in[foot * 2 + i];
      }
    }
  }
  //pthread_mutex_unlock(&solve_mutex);
}

/**
 * Creates a new state error weighting matrix.
 *
 *  @param newQ
 *      the new matrix values for Q in an array.
 */
  template<typename T>
void JumpNLP<T>::SetStateWeights(double * Q_in)
{
  //pthread_mutex_lock(&solve_mutex);
  for (int i = 0; i < NUM_STATES; i++) {
    Q[i] = Q_in[i];
  }
  //pthread_mutex_unlock(&solve_mutex);
}

/**
 * Creates a new input regularization weighting matrix.
 *
 *  @param newR
 *      the new matrix values for R in an array.
 */
  template<typename T>
void JumpNLP<T>::SetInputRegularization(double * R_in)
{
  //pthread_mutex_lock(&solve_mutex);
  int k = 0;
  for (int foot = 0; foot < NUM_FEET; foot++) {
    for (int i = 0; i < 4; i++) {
      R[k] = R_in[i];
      k++;
    }
  }
  //pthread_mutex_unlock(&solve_mutex);
}
template<typename T>
int JumpNLP<T>::GetSolved() {
    //pthread_mutex_lock(&solve_mutex);
    return SOLVED;
    //pthread_mutex_unlock(&solve_mutex);
}
template<typename T>
void JumpNLP<T>::SetRobotParameters(double mass_in, double Inertia_in) {
    //pthread_mutex_lock(&solve_mutex);
    _mass = mass_in;
    _I_yy = Inertia_in;
    //pthread_mutex_unlock(&solve_mutex);
    //cout << "MASS: " << mass << "   |   INERTIA: " << I_body_diagonal[0]<< I_body_diagonal[1]<< I_body_diagonal[2]<< endl;
}


template<typename T>
void JumpNLP<T>::SetEnvironmentParameters(double gravity_in, double mu_in, double * z_g_in) {
  //pthread_mutex_lock(&solve_mutex);

  _gravity = gravity_in;
    _mu = mu_in;
  for (int foot = 0; foot < NUM_FEET; foot++) {
    z_g[foot] = z_g_in[foot];
  }
  //pthread_mutex_unlock(&solve_mutex);
}

template<typename T>
void JumpNLP<T>::SetCurrentState(double * current_state_in) {
  //pthread_mutex_lock(&solve_mutex);
  for (int i = 0; i < NUM_STATES; i++) {
    current_state[i] = current_state_in[i];
  }
  //pthread_mutex_unlock(&solve_mutex);
}

// Get the final solution
template<typename T>
double JumpNLP<T>::GetSolution(int index) {
  //cout << X_final[index-1] << endl;
  return X_final[index - 1];
}
/**
 * Predicted constrints over the given gait-generalized prediction horizon.
 */
template<typename T>
void JumpNLP<T>::Constraints(const double * X, double * constraints) {
  //cout << "NEW ITER" << endl;
  // Run the iterative constraints
  for (int k = 0; k < NUM_PREDICTIONS; k++) {

    // Start indices
    i_x = k * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_f = k * NUM_FEET;
    i_c = k * NUM_CONSTRAINTS;
    i_f1 = (k + 1) * NUM_FEET;
    if (k == (NUM_PREDICTIONS - 1)) {
      // Evaluate the constraints at the final timestep
      Jump2DConstraintsFinal(X + i_x, X + i_u, current_state, 
          contact_state_pred + i_f, dt_pred[k], p_foot_0, _mass, _I_yy, _gravity, _mu,
          constraints + i_c);

    } else {
      i_x1 = (k + 1) * NUM_DECISION_VARS;
      i_u1 = i_x1 + NUM_STATES;

      // Evaluate the constraints at the timestep
      Jump2DConstraints(X + i_x, X + i_u, X + i_x1, contact_state_pred+i_f, 
          dt_pred[k], p_foot_0, _mass, _I_yy, _gravity, _mu, 
          constraints + i_c);
    }
    /*cout << X[i_x+2];
      for (int i = 0; i < NUM_FEET; i++) {
      cout << " ground: " << z_g[i] << " r's: " << X[i_u+2 + i*6] << " " ;
      }
      cout << endl;*/
    /*
       for (int i = 0; i < NUM_CONSTRAINTS; i++) {
       cout << constraints[i] << " " ;
       }
       cout << "\n\n" << endl;
       */
  }
}


/**
 * Predicted lagrangian hessian over the given gait-generalized prediction horizon.
 */
template<typename T>
void JumpNLP<T>::LagrangianHessian(const double * X, double * hessian, Number obj_factor, const Number * lambda) {

  (void)X;
  (void)lambda;

  // Incremental cost over the predicted gait horizon
  for (int k = 0; k < NUM_PREDICTIONS; k++) {

    // Start indices
    i_x = k * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_xd = k * NUM_STATES;
    i_c = k * NUM_CONSTRAINTS;
    i_h = k * NNZ_H;
    i_f = k * NUM_FEET;
    i_Q = k * NUM_STATES;
    i_R = k * NUM_INPUTS;

    if (k == (NUM_PREDICTIONS - 1)) {
      // Calculate the iterative Hessian of the Lagrangian at the lfinal time
      Jump2DLagrangianHessianFinal(Q + i_Q, R + i_R, obj_factor, hessian + i_h);
    } else {
      // Calculate the iterative Hessian of the Lagrangian over the prediction horizon
      Jump2DLagrangianHessian(Q + i_Q, R + i_R, obj_factor, hessian + i_h);
    }
  }
}


/**
 * Predicted cost over the given gait-generalized prediction horizon.
 */
template<typename T>
void JumpNLP<T>::ConstraintJacobian(const double * X, double * constraint_jacobian) {
  (void)X;
  // Incremental cost over the predicted gait horizon
  for (int k = 0; k < NUM_PREDICTIONS; k++) {
    i_jc = k * NNZ_J_G;
    i_x = k * NUM_DECISION_VARS;
    i_u = i_x + NUM_STATES;
    i_f = k * NUM_FEET;
    i_c = k * NUM_CONSTRAINTS;
    if (k == (NUM_PREDICTIONS - 1)) {
      i_x1 = i_x;
      i_u1 = i_u;
    } else {
      i_x1 = (k + 1) * NUM_DECISION_VARS;
      i_u1 = i_x1 + NUM_STATES;
    }
    i_f1 = (k + 1) * NUM_FEET;

    if (k == (NUM_PREDICTIONS - 1)) {
      // Final constraints only dependent on contact state and friction
      Jump2DConstraintJacobianFinal(contact_state_pred + i_f, dt_pred[k], p_foot_0, 
          _mass, _I_yy, _mu,
          constraint_jacobian + i_jc);
    } else {
      Jump2DConstraintJacobian(
          contact_state_pred + i_f, dt_pred[k], p_foot_0,
          _mass, _I_yy, _mu,
          constraint_jacobian + i_jc);

    }
  }
}


template class JumpNLP<double>;
template class JumpNLP<float>;
#endif
