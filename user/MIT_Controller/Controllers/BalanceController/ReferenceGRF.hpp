#ifndef REFERENCEGRF_H
#define REFERENCEGRF_H
/*
+References:
+   [R1] M. Focchi, A. del Prete, I. Havoutis, R. Featherstone, D. G. Caldwell, and C. Semini. High-slope terrain
+   locomotion for torque-controlled quadruped robots. Autonomous Robots, 2016.
+   [R2] R. M. Murray, S. S. Sastry, and L. Zexiang. A Mathematical Introduction to Robotic Manipulation. CRC
+   Press, Inc., Boca Raton, FL, USA, 1st edition, 1994.
+Cheetah-3-Documentation-Control:
+   [C1] balanceController.pdf
+   qpOASES variables are terminated with _qpOASES
+*/
   
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <eigen3/Eigen/Dense>
#include <qpOASES.hpp>
#include <lcm/lcm-cpp.hpp>
#include "sim_command_t.hpp"
#include "qp_controller_data_t.hpp"


static const int NUM_VARIABLES_QP_DES = 4;
static const int NUM_CONSTRAINTS_QP_DES = 4;
static const int NUM_CONTACT_POINTS_DES = 4;
static const int NUM_VARIABLES_PER_FOOT_DES = 1;
static const int NUM_CONSTRAINTS_PER_FOOT_DES = 1;

//static const double PI_CONST = 3.1415;
static const double NEGATIVE_NUMBER_DES = -1000000.0;
static const double POSITIVE_NUMBER_DES =  1000000.0;

using namespace Eigen;
using namespace qpOASES;

class ReferenceGRF
{
   public:
      ReferenceGRF();
      ~ReferenceGRF(){};

      // use new kinematics measurements to update QP
      void updateProblemData(double* p_feet_desired_in,double* p_des); // new

      // Min/max forces to set QP constraints
      void SetContactData(double* contact_state_in,double* min_forces_in,double* max_forces_in, double threshold, int stance_legs_in);

      // calculate the QP, return solution
      void solveQP_nonThreaded(double* xOpt);

      // configure gains, QP weights, force limits, world parameters
      void set_RobotLimits();
      void set_worldData();
      void set_QPWeights();
      void set_mass(double mass_in);
      void set_alpha_control(double alpha_control_in);


   private:
      /* Fixed-Size qpOASES data */           
      QProblemB QProblemObj_qpOASES;  
      int_t nWSR_qpOASES = 100;
      int_t nWSR_fixed = 100;
      real_t cpu_time;
      real_t cpu_time_fixed;
      int_t qp_exit_flag;
      int nWSR_initial;
      double cpu_time_initial;
      double xOpt_local[4];
      double qp_not_init;

      Bounds guessedBounds;
      Constraints guessedConstraints;
      
      /* QP variables for HJB optimization */
      real_t    H_qpOASES[NUM_VARIABLES_QP_DES*NUM_VARIABLES_QP_DES];
      real_t    g_qpOASES[NUM_VARIABLES_QP_DES];
      real_t   lb_qpOASES[NUM_VARIABLES_QP_DES];
      real_t   ub_qpOASES[NUM_VARIABLES_QP_DES];  
      real_t xOpt_qpOASES[NUM_VARIABLES_QP_DES];
      real_t yOpt_qpOASES[NUM_VARIABLES_QP_DES+NUM_CONSTRAINTS_QP_DES];   
      real_t xOpt_initialGuess[NUM_VARIABLES_QP_DES];

      /* Eigen Variables that Match qpOASES variables */
      Eigen::MatrixXd    H_eigen;
      Eigen::MatrixXd    g_eigen;
      Eigen::VectorXd xOpt_eigen;
      Eigen::VectorXd yOpt_eigen;

      /* Robot control variables used to construct QP matrices, see (5) and (6) of [R1] */
      Eigen::MatrixXd A_control;            
      Eigen::VectorXd b_control;
      double alpha_control;

      /* Model and World parameters and force limits */
      double mass;
      Eigen::VectorXd gravity;

      Eigen::VectorXd minNormalForces_feet;
      Eigen::VectorXd maxNormalForces_feet;

      /* Foot Contact Information, 1 is on the ground,  */
      Eigen::VectorXd contact_state;

      /* Desired Kinematics */
      Eigen::VectorXd x_COM_world_desired;
      Eigen::MatrixXd p_feet_desired; //new 

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
      void calc_H_qpOASES();
      void calc_g_qpOASES();     
      void calc_lb_ub_qpOASES();

      /* Utility Functions */
      void copy_Eigen_to_real_t(real_t* target, Eigen::MatrixXd &source, int nRows, int nCols );                  
      void copy_Array_to_Eigen(Eigen::VectorXd &target, double* source, int len, int startIndex);
      void copy_Array_to_Eigen(Eigen::MatrixXd &target, double* source, int len, int startIndex);
      void copy_real_t_to_Eigen(Eigen::VectorXd &target, real_t* source, int len);
      void crossMatrix(Eigen::MatrixXd &R, const Eigen::VectorXd &omega);
};

#endif
