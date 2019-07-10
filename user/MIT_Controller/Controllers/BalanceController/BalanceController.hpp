#ifndef BALANCECONTROLLER_H
#define BALANCECONTROLLER_H
/*

References:
   [R1] M. Focchi, A. del Prete, I. Havoutis, R. Featherstone, D. G. Caldwell,
and C. Semini. High-slope terrain locomotion for torque-controlled quadruped
robots. Autonomous Robots, 2016.

   [R2] R. M. Murray, S. S. Sastry, and L. Zexiang. A Mathematical Introduction
to Robotic Manipulation. CRC Press, Inc., Boca Raton, FL, USA, 1st edition,
1994.

Cheetah-3-Documentation-Control:
   [C1] balanceController.pdf

   qpOASES variables are terminated with _qpOASES
*/

#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <eigen3/Eigen/Dense>
#include <lcm/lcm-cpp.hpp>
#include <qpOASES.hpp>
#include "qp_controller_data_t.hpp"
#include "sim_command_t.hpp"

static const int NUM_VARIABLES_QP = 12;
static const int NUM_CONSTRAINTS_QP = 20;
static const int NUM_CONTACT_POINTS = 4;
static const int NUM_VARIABLES_PER_FOOT = 3;
static const int NUM_CONSTRAINTS_PER_FOOT = 5;

// static const double PI_CONST = 3.1415;
static const double NEGATIVE_NUMBER = -1000000.0;
static const double POSITIVE_NUMBER = 1000000.0;

using namespace Eigen;
using namespace qpOASES;

class BalanceController {
 public:
  BalanceController();
  ~BalanceController(){};

  void testFunction();

  // use new kinematics measurements to update QP
  void updateProblemData(double* xfb_in, double* p_feet_in, double* p_des,
                         double* p_act, double* v_des, double* v_act,
                         double* O_err, double yaw_act_in);

  void SetContactData(double* contact_state_in, double* min_forces_in,
                      double* max_forces_in);

  // calculate the QP, return solution
  void solveQP(double* xOpt);
  void solveQP_nonThreaded(double* xOpt);

  // update desired COM and orientation setpoints
  void set_desiredTrajectoryData(double* rpy_des_in, double* p_des_in,
                                 double* omegab_des_in, double* v_des_in);
  void set_PDgains(double* Kp_COM_in, double* Kd_COM_in, double* Kp_Base_in,
                   double* Kd_Base_in);
  void set_QP_options(double use_hard_constraint_pitch_in);

  // configure gains, QP weights, force limits, world parameters
  void set_RobotLimits();
  void set_worldData();
  void set_PDgains();
  void set_wrench_weights(double* COM_weights_in, double* Base_weights_in);
  void set_QPWeights();
  void set_friction(double mu_in);
  void set_alpha_control(double alpha_control_in);
  void set_mass(double mass_in);
  void set_desired_swing_pos(double* pFeet_des_in);
  void set_actual_swing_pos(double* pFeet_act_in);

  // print QP matrices (real_t)
  void print_QPData();

  void verifyModel(double* vbd_command);
  void set_base_support_flag(double sflag);
  void publish_data_lcm();

 private:
  lcm::LCM* lcm;
  qp_controller_data_t qp_controller_data, qp_controller_data_publish;
  sim_command_t command;

  /* Fixed-Size qpOASES data */
  QProblem QProblemObj_qpOASES;

  int_t nWSR_qpOASES = 100;
  int_t nWSR_fixed = 100;

  real_t cpu_time;
  real_t cpu_time_fixed;

  int_t qp_exit_flag;

  int nWSR_initial;
  double cpu_time_initial;

  double xOpt_local[12];
  double qp_not_init;

  Bounds guessedBounds;
  Constraints guessedConstraints;

  real_t H_qpOASES[NUM_VARIABLES_QP * NUM_VARIABLES_QP];
  real_t A_qpOASES[NUM_CONSTRAINTS_QP * NUM_VARIABLES_QP];
  real_t g_qpOASES[NUM_VARIABLES_QP];
  real_t lb_qpOASES[NUM_VARIABLES_QP];
  real_t ub_qpOASES[NUM_VARIABLES_QP];
  real_t lbA_qpOASES[NUM_CONSTRAINTS_QP];
  real_t ubA_qpOASES[NUM_CONSTRAINTS_QP];
  real_t xOpt_qpOASES[NUM_VARIABLES_QP];
  real_t yOpt_qpOASES[NUM_VARIABLES_QP + NUM_CONSTRAINTS_QP];

  real_t xOpt_initialGuess[NUM_VARIABLES_QP];

  /* Eigen Variables that Match qpOASES variables */
  Eigen::MatrixXd H_eigen;
  Eigen::MatrixXd A_eigen;
  Eigen::MatrixXd g_eigen;
  Eigen::VectorXd xOpt_eigen;
  Eigen::VectorXd yOpt_eigen;

  /* Robot control variables used to construct QP matrices, see (5) and (6) of
   * [R1] */
  Eigen::MatrixXd A_control;
  Eigen::MatrixXd S_control;
  Eigen::MatrixXd W_control;
  Eigen::MatrixXd C_control;
  Eigen::VectorXd b_control;
  Eigen::VectorXd b_control_Opt;

  Eigen::VectorXd C_times_f_control;

  double alpha_control;

  /* Centroidal control PD gains and variables, see (3) and (4) of [R1] */
  double Kp_COMx, Kp_COMy, Kp_COMz;
  double Kd_COMx, Kd_COMy, Kd_COMz;

  double Kp_Base_roll, Kp_Base_pitch, Kp_Base_yaw;
  double Kd_Base_roll, Kd_Base_pitch, Kd_Base_yaw;

  double use_hard_constraint_pitch;

  /* Model and World parameters and force limits */
  double mass;
  double inertia;
  Eigen::MatrixXd Ig;

  double mu_friction;
  Eigen::VectorXd gravity;

  Eigen::VectorXd minNormalForces_feet;
  Eigen::VectorXd maxNormalForces_feet;

  Eigen::VectorXd direction_normal_flatGround;
  Eigen::VectorXd direction_tangential_flatGround;

  /* Foot Contact Information, 1 is on the ground,  */
  Eigen::VectorXd contact_state;

  double yaw_act;

  /* Actual Kinematics*/
  Eigen::VectorXd x_COM_world;
  Eigen::VectorXd xdot_COM_world;
  Eigen::VectorXd omega_b_world;
  Eigen::VectorXd quat_b_world;
  Eigen::MatrixXd R_b_world;
  Eigen::MatrixXd p_feet;

  Eigen::MatrixXd R_yaw_act;

  /* Desired Kinematics */
  Eigen::VectorXd x_COM_world_desired;
  Eigen::VectorXd xdot_COM_world_desired;
  Eigen::VectorXd xddot_COM_world_desired;
  Eigen::VectorXd omega_b_world_desired;
  Eigen::VectorXd omegadot_b_world_desired;
  Eigen::MatrixXd R_b_world_desired;

  // Error coordinates
  Eigen::VectorXd error_x_rotated;
  Eigen::VectorXd error_dx_rotated;
  Eigen::VectorXd error_theta_rotated;
  Eigen::VectorXd error_dtheta_rotated;

  Eigen::VectorXd vbd_command_eigen;

  Eigen::VectorXd orientation_error;

  /* Temporary, Internal Matrices */
  Eigen::MatrixXd omegaHat;
  Eigen::MatrixXd tempSkewMatrix3;
  Eigen::VectorXd tempVector3;

  /* Interal QP management data */
  bool QPFinished;

  Eigen::VectorXd xOptPrev;
  Eigen::VectorXd yOptPrev;

  /* Interface Functions */
  bool getQPFinished();

  void update_A_control();
  void calc_PDcontrol();

  void calc_H_qpOASES();
  void calc_A_qpOASES();
  void calc_g_qpOASES();
  void calc_lb_ub_qpOASES();
  void calc_lbA_ubA_qpOASES();
  void update_log_variables(double* p_des, double* p_act, double* v_des,
                            double* v_act, double* O_err);
  void calc_constraint_check();

  /* Utility Functions */
  void copy_Eigen_to_real_t(real_t* target, Eigen::MatrixXd& source, int nRows,
                            int nCols);
  void copy_Eigen_to_double(double* target, Eigen::VectorXd& source,
                            int length);
  void copy_Array_to_Eigen(Eigen::VectorXd& target, double* source, int len,
                           int startIndex);
  void copy_Array_to_Eigen(Eigen::MatrixXd& target, double* source, int len,
                           int startIndex);
  void copy_real_t_to_Eigen(Eigen::VectorXd& target, real_t* source, int len);

  void print_real_t(real_t* matrix, int nRows, int nCols);

  void matrixExpOmegaCross(const Eigen::VectorXd& omega, Eigen::MatrixXd& R);
  void matrixLogRot(const Eigen::MatrixXd& R, Eigen::VectorXd& omega);
  void crossMatrix(Eigen::MatrixXd& R, const Eigen::VectorXd& omega);
  void quaternion_to_rotationMatrix(Eigen::MatrixXd& R, Eigen::VectorXd& quat);

  void rpyToR(Eigen::MatrixXd& R, double* rpy_in);
};

#endif
