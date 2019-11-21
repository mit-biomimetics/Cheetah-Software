/*

References:
   [R1] M. Focchi, A. del Prete, I. Havoutis, R. Featherstone, D. G. Caldwell, and C. Semini. High-slope terrain
   locomotion for torque-controlled quadruped robots. Autonomous Robots, 2016.

   [R2] R. M. Murray, S. S. Sastry, and L. Zexiang. A Mathematical Introduction to Robotic Manipulation. CRC
   Press, Inc., Boca Raton, FL, USA, 1st edition, 1994.

Cheetah-3-Documentation-Control:
   [C1] BalanceControllerVBL.pdf

   qpOASES variables are terminated with _qpOASES
*/


#include "ReferenceGRF.hpp"
#include<iostream>
#include<iomanip>
#include<cstdlib>

// #include "Sim_Utils.h"

ReferenceGRF::ReferenceGRF() :  QProblemObj_qpOASES(NUM_VARIABLES_QP_DES)
{
   // *!! It seems like this next block of 5 lines does not actually do anything, 
   // I had to set print level in the OptimizationThread class
   Options options;
   options.printLevel = PL_NONE;
   QProblemObj_qpOASES.setOptions( options );    
   QProblemObj_qpOASES.setPrintLevel(PL_NONE);
   QPFinished = false;


   // Eigen QP matrices
   H_eigen.resize(NUM_VARIABLES_QP_DES,NUM_VARIABLES_QP_DES);
   g_eigen.resize(NUM_VARIABLES_QP_DES,1);
   xOpt_eigen.resize(NUM_VARIABLES_QP_DES,1);
   yOpt_eigen.resize(NUM_VARIABLES_QP_DES+NUM_CONSTRAINTS_QP_DES,1);

   /* Model and World parameters and force limits */
   mass = 41;
   gravity.resize(3,1);

   /* Desired Kinematics */
   x_COM_world_desired.resize(3,1);
   p_feet_desired.resize(3,4);

   /* Initialize to all feet on the ground */
   contact_state.resize(4,1);
   contact_state << 1,1,1,1;

   minNormalForces_feet.resize(4,1);
   maxNormalForces_feet.resize(4,1);

   /* Temporary, Internal Matrices */
   omegaHat.resize(3,3);   
   tempSkewMatrix3.resize(3,3);
   tempVector3.resize(3,1);

   tempSkewMatrix3.setZero();
   tempVector3.setZero();
   
   /* Robot Control variables for QP Matrices */
   A_control.resize(3,NUM_VARIABLES_QP_DES);
   b_control.resize(3,1);

   xOptPrev.setZero(NUM_VARIABLES_QP_DES);
   yOptPrev.setZero(NUM_VARIABLES_QP_DES+NUM_CONSTRAINTS_QP_DES);

   /* Initialize Optimization variables */
   for(int i = 0; i < NUM_VARIABLES_QP_DES; i++)
   {
      xOpt_qpOASES[i] = 0.0;
      xOpt_initialGuess[i] = 99.0;
   }

   for(int i = 0; i < NUM_VARIABLES_QP_DES+NUM_CONSTRAINTS_QP_DES; i++)
   {
      yOpt_qpOASES[i] = 0.0;
   }

   set_RobotLimits();
   set_worldData();
   set_QPWeights();      

   cpu_time = 0.001;
   cpu_time_fixed = 0.001;
   qp_exit_flag = -1.0;

   qp_not_init = 1.0;
}



/* --------------- Primary Interface ------------------- */

void ReferenceGRF::updateProblemData(double* p_feet_desired_in,double* p_des_in)                           
{   

  // Unpack inputs
  copy_Array_to_Eigen(x_COM_world_desired, p_des_in, 3,0);
  copy_Array_to_Eigen(p_feet_desired, p_feet_desired_in, 12,0);

   // Update Centroidal Robot dynamics
   update_A_control(); // set A control and b control based on dyanmics

   // Compute QP Problem data
   calc_H_qpOASES();
   calc_g_qpOASES(); 

   cpu_time = cpu_time_fixed;
   nWSR_qpOASES = nWSR_fixed;  

}

void ReferenceGRF::SetContactData( double* contact_state_in,
                           double* min_forces_in,
                           double* max_forces_in,
                           double threshold,
                           int stance_legs_in)                                                  
{   
   for(int leg = 0; leg < 4; leg++){
    if(contact_state_in[leg] > threshold)
      contact_state(leg) = 1.0;
    else
      contact_state(leg) = 0.0;
   }

   if(stance_legs_in == 5){
    contact_state(0) = 0.0;
    contact_state(1) = 1.0;
    contact_state(2) = 1.0;
    contact_state(3) = 0.0;
   }

   // Unpack inputs
   copy_Array_to_Eigen(minNormalForces_feet, min_forces_in,4,0);
   copy_Array_to_Eigen(maxNormalForces_feet, max_forces_in,4,0);
   
   calc_lb_ub_qpOASES();
}


void ReferenceGRF::solveQP_nonThreaded(double* xOpt)
{
  // &cpu_time
  if(qp_not_init==1.0)
  {
    qp_exit_flag = QProblemObj_qpOASES.init(H_qpOASES,g_qpOASES,lb_qpOASES,ub_qpOASES,nWSR_qpOASES, &cpu_time, xOpt_initialGuess); 
    qp_not_init = 0.0;

    nWSR_initial = nWSR_qpOASES;
    cpu_time_initial = cpu_time;
  }else
  {
    qp_exit_flag = QProblemObj_qpOASES.init(H_qpOASES,g_qpOASES,lb_qpOASES,ub_qpOASES,nWSR_qpOASES,&cpu_time,xOpt_qpOASES);
  }
   
  QProblemObj_qpOASES.getPrimalSolution(xOpt_qpOASES);
  // QProblemObj_qpOASES.getDualSolution( yOpt_qpOASES );

  // QProblemObj_qpOASES.getBounds( guessedBounds );
  // QProblemObj_qpOASES.getConstraints( guessedConstraints );

  copy_real_t_to_Eigen(xOpt_eigen, xOpt_qpOASES, 4);

  xOptPrev = xOpt_eigen;

   
   for(int i = 0; i < NUM_VARIABLES_QP_DES; i++)
   {
      xOpt[i] = xOpt_eigen(i);
   }

  QProblemObj_qpOASES.reset();
}

/* --------------- QP Matrices and Problem Data ------------ */
void ReferenceGRF::update_A_control()
{

   // Update the A matrix and b vector in the controller notation A*f = b
  A_control.block<1,4>(0,0) << 1,1,1,1;
  A_control.block<1,4>(1,0) << p_feet_desired.row(1);
  A_control.block<1,4>(2,0) << -1*p_feet_desired.row(0);
  b_control << mass*gravity(2),0,0;
}


void ReferenceGRF::calc_H_qpOASES()
{
  // Use LQR matrices to compute the QP cost matrix H
  Eigen::MatrixXd IDNTY;IDNTY.setIdentity(NUM_VARIABLES_QP_DES,NUM_VARIABLES_QP_DES);
  H_eigen = 2*(A_control.transpose()*A_control+alpha_control*IDNTY);

  // Copy to real_t array (qpOASES data type)
   copy_Eigen_to_real_t(H_qpOASES, H_eigen, NUM_VARIABLES_QP_DES, NUM_VARIABLES_QP_DES);
}

void ReferenceGRF::calc_g_qpOASES()
{  
   // Use LQR matrices to compute the QP cost vector g
   g_eigen = -2*A_control.transpose()*b_control;

   // Copy to real_t array (qpOASES data type)
   copy_Eigen_to_real_t(g_qpOASES, g_eigen, NUM_VARIABLES_QP_DES, 1);
}


void ReferenceGRF::calc_lb_ub_qpOASES()
{
   for(int i = 0; i < NUM_VARIABLES_QP_DES; i++)
   {
      lb_qpOASES[i] = contact_state(i)*NEGATIVE_NUMBER_DES;
      ub_qpOASES[i] = contact_state(i)*POSITIVE_NUMBER_DES;    
   }

}


/* ------------ Set Parameter Values -------------- */
void ReferenceGRF::set_worldData()
{
   gravity << 0,0,9.81;
}

void ReferenceGRF::set_mass(double mass_in)
{
   mass = mass_in;
}

void ReferenceGRF::set_QPWeights()
{
   alpha_control = .01;
}

void ReferenceGRF::set_RobotLimits()
{
   minNormalForces_feet << 10, 10, 10, 10;
   maxNormalForces_feet << 160, 160, 160, 160;
}

void ReferenceGRF::set_alpha_control(double alpha_control_in)
{
   alpha_control = alpha_control_in;
}

/* ------------ Utilities -------------- */


bool ReferenceGRF::getQPFinished()
{
   return QPFinished;
}

void ReferenceGRF::copy_Eigen_to_real_t(real_t* target, Eigen::MatrixXd &source, int nRows, int nCols )
{
   int count = 0;

   // Strange Behavior: Eigen matrix matrix(count) is stored by columns (not rows) 
   for(int i = 0; i < nRows; i++)
   {      
      for (int j = 0; j < nCols; j++)
      {
         target[count] = source(i,j);
         count++;
      }
   }  
}


void ReferenceGRF::copy_Array_to_Eigen(Eigen::VectorXd &target, double* source, int len, int startIndex)
{
   for(int i = 0; i < len; i++)
   {
      target(i) = source[i+startIndex];
   }
}

void ReferenceGRF::copy_Array_to_Eigen(Eigen::MatrixXd &target, double* source, int len, int startIndex)
{
   for(int i = 0; i < len; i++)
   {
      target(i) = source[i+startIndex];
   }
}

void ReferenceGRF::copy_real_t_to_Eigen(Eigen::VectorXd &target, real_t* source, int len)
{
   for(int i = 0; i < len; i++)
   {
      target(i) = source[i];
   }
}

void ReferenceGRF::crossMatrix(Eigen::MatrixXd &R, const Eigen::VectorXd &omega)
{
   R(0,1) = -omega(2);
   R(0,2) = omega(1);
   R(1,0) = omega(2);
   R(1,2) = -omega(0);
   R(2,0) = -omega(1);
   R(2,1) = omega(0);
}
