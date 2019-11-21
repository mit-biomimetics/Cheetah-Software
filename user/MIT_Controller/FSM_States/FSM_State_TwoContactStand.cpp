/*============================= Two Contact Stand ==============================*/
/**
 * Transitionary state that is called for the robot to stand up into
 * baland control mode.
 */

#include "FSM_State_TwoContactStand.h"
#include <Utilities/Utilities_print.h>
#include <iostream>
#include <iomanip>
#include <fstream>

//#define DRAW_DEBUG


/**
 * Constructor for the FSM State that passes in state specific info to
 * the generic FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_TwoContactStand<T>::FSM_State_TwoContactStand(ControlFSMData<T>* _controlFSMData)
    : FSM_State<T>(_controlFSMData, FSM_StateName::TWO_CONTACT_STAND, "TWO_CONTACT_STAND"){
  // Post control safety checks
  this->checkPDesFoot = false;
  this->checkForceFeedForward = false;

  // Build Quadruped Model
  model = _controlFSMData->_quadruped->buildModel();
  _wbc_ctrl = new LocomotionCtrl<T>(_controlFSMData->_quadruped->buildModel());
  _wbc_data = new LocomotionCtrlData<T>();

}

template <typename T>
void FSM_State_TwoContactStand<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();

  // Reset iteration counter
  iter = 0;

  // Reset E-stop
  ESTOP = false;

  // Remove impedance control for all joints
  impedance_kp << 0.0, 0.0, 0.0;
  impedance_kd << 0.0, 0.0, 0.0;

  // Initial Contact Phase
  conPhase << 0.5, 0.5, 0.5, 0.5;

}

/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_TwoContactStand<T>::run() {
  // Count iterations
  iter++;

  // Get current state
  get_current_state();

  // Get the foot locations relative to COM
  get_foot_locations();

  // Get desired state from gamepad controls
  get_desired_state();

  // Set Control Parameters
  set_ctrl_params();

  // Compute Reference Control Inputs - only if desired xy pos has changed
  if (p_des_prev[0] != p_des[0] || p_des_prev[1] != p_des[1]) {
    refGRF.SetContactData(contactStateScheduled, minForces, maxForces, threshold, this->_data->userParameters->stance_legs);
    refGRF.updateProblemData(pFeet_yaw_des, p_des);
    refGRF.solveQP_nonThreaded(f_ref_z);
    for (int leg = 0; leg < 4; leg++)
      f_ref_world[3*leg+2] = f_ref_z[leg];
    balanceControllerVBL.set_reference_GRF(f_ref_world);
  }

  // Solve Balance control QP
  balanceControllerVBL.set_desiredTrajectoryData(rpy, p_des, omegaDes, v_des);
  balanceControllerVBL.SetContactData(contactStateScheduled, minForces, maxForces, threshold, this->_data->userParameters->stance_legs);
  balanceControllerVBL.updateProblemData(se_xfb, pFeet_yaw, pFeet_yaw_des, rpy, rpy_yaw_offset);
  balanceControllerVBL.solveQP_nonThreaded(fOpt);
  balanceControllerVBL.publish_data_lcm();

  // Compute feedforward forces
  compute_force_ff();

  // Send commands to leg controller or wbc
  if(this->_data->userParameters->use_wbc > 0.9 && ESTOP == false){
    _wbc_data->pBody_des = pBody_des;
    _wbc_data->vBody_des = vBody_des;
    _wbc_data->aBody_des = aBody_des;

    _wbc_data->pBody_RPY_des = pBody_RPY_des;
    _wbc_data->vBody_Ori_des = vBody_Ori_des;
    
    for(size_t i(0); i<4; ++i){
      _wbc_data->pFoot_des[i] = pFoot_des[i];
      _wbc_data->vFoot_des[i] = vFoot_des[i];
      _wbc_data->aFoot_des[i] = aFoot_des[i];
      _wbc_data->Fr_des[i] = Fr_des[i];
    }
    _wbc_data->contact_state = contact_state;
    _wbc_ctrl->run(_wbc_data, *this->_data);

  } else {
    // Zero the leg controller
    this->_data->_legController->commands->zero();

    for (int leg = 0; leg < 4; leg++) {
      if(contactStateScheduled[leg] > 0.0){
        // Impedance Control
        this->cartesianImpedanceControl(leg, this->footstepLocations.col(leg), Vec3<T>::Zero(),impedance_kp,impedance_kd);
      
        // Force and Joint Torque control
        this->_data->_legController->commands[leg].forceFeedForward = this->footFeedForwardForces.col(leg);
        this->_data->_legController->commands[leg].tauFeedForward = this->jointFeedForwardTorques.col(leg);
        
      } else {
        this->liftLeg(leg, q_lift_leg, qd_lift_leg);

      }
    }
  }

  // Visualize for debugging
  #ifdef DRAW_DEBUG
    debug_visualization();
  #endif

  // Set the contact phase for the state estimator
  this->_data->_stateEstimator->setContactPhase(conPhase);

  // Update previous desired posiion
  for (int i = 0; i < 2; i++)
    p_des_prev[i] = p_des[i];

}

/**
* Computes the feedforces using the output from the optimization
**/
template <typename T>
void FSM_State_TwoContactStand<T>::compute_force_ff() {
  
  // Check for emergency stop if orientation error too large
  for(int i = 0; i < 3; i++){
    if(fabs(rpy[i]-rpy_yaw_offset[i])>0.85)
      ESTOP = true;
  }

  // Feed forward forces for legs in contact with the ground & PD control for legs not in contact
  qd_lift_leg << 0.0, 0.0, 0.0;
  if(ESTOP){
    this->jointFeedForwardTorques = Mat34<float>::Zero();   // feed forward joint torques
    this->footFeedForwardForces = Mat34<float>::Zero();     // feedforward forces at the feet
    for(int i = 0; i < 4; i++)
      Fr_des[i] = Vec3<float>::Zero();
  } else{
  for (int leg = 0; leg < 4; leg++) {
      this->footFeedForwardForces.col(leg) << (T)fOpt[leg * 3], (T)fOpt[leg * 3 + 1], (T)fOpt[leg * 3 + 2];
      this->jointFeedForwardTorques.col(leg) << (T)G_ff[6+3*leg], (T)G_ff[6+3*leg+1], G_ff[6+3*leg+2];
    }
  }

  // For WBC
  for(int i = 0; i < 3; i++){
    pBody_des[i] = p_des[i];
    vBody_des[i] = v_des[i];
    aBody_des[i] = 0.0;
    pBody_RPY_des[i] = rpy[i];
    vBody_Ori_des[i] = omegaDes[i];
  }
  pBody_RPY_des[2] = rpy[2] + ini_yaw; // remove initial yaw compensation

  // need to compute Jleg
  Vec3<float> force_ff;
  Vec3<float> grav_comp;
  for(int leg = 0; leg < 4; leg++){
    force_ff = {(float) -fOpt[3*leg],(float) -fOpt[3*leg+1],(float) -fOpt[3*leg+2]};
    grav_comp = {(float) G_ff[6+3*leg],(float) G_ff[6+3*leg+1],(float) G_ff[6+3*leg+2]};
    computeLegJacobianAndPosition(**&this->_data->_quadruped, this->_data->_legController->datas[leg].q, &Jleg, (Vec3<T>*)nullptr, leg);
    Fr_des[leg] = this->_data->_stateEstimator->getResult().rBody.transpose() * (force_ff + Jleg.transpose().inverse()*grav_comp);
    pFoot_des[leg] = {(float) pFeet_world_des[3*leg],(float) pFeet_world_des[3*leg+1],(float) pFeet_world_des[3*leg+2]};
    vFoot_des[leg] = {0.0, 0.0, 0.0};
    aFoot_des[leg] = {0.0, 0.0, 0.0};
    contact_state[leg] = conPhase[leg]; 
  }


}


/**
* Get's the current state of the robot from state estimator. Accounts for non-zero initial yaw
* Also finds the position of the CoM relative to the body frame
**/
template <typename T>
void FSM_State_TwoContactStand<T>::get_current_state() {
  
  // Get orientation from state estimator
  for (int i = 0; i < 4; i++) {
    quat_act[i] = (double)this->_data->_stateEstimator->getResult().orientation(i);
  }

  // Get current state from state estimator
  for (int i = 0; i < 3; i++) {
    p_act[i] = (double)this->_data->_stateEstimator->getResult().position(i);
    se_xfb[7 + i] = (double)this->_data->_stateEstimator->getResult().omegaBody(i);
    se_xfb[10 + i] = (double)this->_data->_stateEstimator->getResult().vBody(i);
  }

  // Convert quaternions to RPY that accounts for initial yaw offset
  quatToEuler(quat_act, rpy_yaw_offset);

  // Account for initial yaw
  if (iter < 5)
    ini_yaw = rpy_yaw_offset[2];
  rpy_yaw_offset[2] -= ini_yaw;

  // Convert back to quaternions w/ initial yaw accounted for
  eulerToQuat(rpy_yaw_offset, quat_yaw_offset);
  for (int i = 0; i < 4; i++) {
    se_xfb[i] = quat_yaw_offset[i];
  }

  // Update state for the quadruped model using SE
  for (int i = 0; i < 3; i++){
    state.bodyOrientation[i] = quat_yaw_offset[i];
    state.bodyPosition[i] = this->_data->_stateEstimator->getResult().position(i);
    state.bodyVelocity[i] = this->_data->_stateEstimator->getResult().omegaBody(i);
    state.bodyVelocity[i+3] = this->_data->_stateEstimator->getResult().vBody(i);
  }
  state.bodyOrientation[3] = quat_yaw_offset[3];
  state.q.setZero(12);
  state.qd.setZero(12);

  for (int i = 0; i < 4; ++i) {
   state.q(3*i+0) = this->_data->_legController->datas[i].q[0];
   state.q(3*i+1) = this->_data->_legController->datas[i].q[1];
   state.q(3*i+2) = this->_data->_legController->datas[i].q[2];
   state.qd(3*i+0)= this->_data->_legController->datas[i].qd[0];
   state.qd(3*i+1)= this->_data->_legController->datas[i].qd[1];
   state.qd(3*i+2)= this->_data->_legController->datas[i].qd[2];
  }
  model.setState(state);
  H = model.massMatrix();
  G_ff = model.generalizedGravityForce();

  // CoM relative to body frame
  mass_in = H(3,3);
  c_body[0] = H(2,4)/mass_in;
  c_body[1] = H(0,5)/mass_in;
  c_body[2] = H(1,3)/mass_in;

  // CoM relative to world frame
  rpyToR(rBody_yaw, rpy_yaw_offset);
  c_world = rBody_yaw * c_body;

  // Position of CoM in world frame
  for (int i = 0; i < 3; i++)
    p_COM[i] = p_act[i]+c_world[i]; // assumes C is expressed in world frame

  // Use the COM pos instead of the body pos for balance control
  for (int i = 0; i < 3; i++) {
    se_xfb[i+4] = p_COM[i]; 
    p_act[i] = p_COM[i];
  }

}

/**
* Set control parameters for the balance controller
*/
template <typename T>
void FSM_State_TwoContactStand<T>::set_ctrl_params() {

  // LQR Weights
  for (int i = 0; i < 3; i++) {
    x_weights[i] = this->_data->userParameters->Q_pos[i];
    xdot_weights[i] = this->_data->userParameters->Q_vel[i];
    R_weights[i] = this->_data->userParameters->Q_ori[i];
    omega_weights[i] = this->_data->userParameters->Q_ang[i];
  }
  alpha_control = this->_data->userParameters->R_control;
  beta_control = this->_data->userParameters->R_prev;
  balanceControllerVBL.set_LQR_weights(x_weights,xdot_weights,R_weights,omega_weights,alpha_control,beta_control);
  refGRF.set_alpha_control(0.01);

  // Friction & inertia
  balanceControllerVBL.set_friction(mu_ctrl);
  balanceControllerVBL.set_mass(mass_in);
  refGRF.set_mass(mass_in);
  if (this->_data->_quadruped->_robotType == RobotType::CHEETAH_3)
    balanceControllerVBL.set_inertia(0.35, 2.10, 2.10);
  else if (this->_data->_quadruped->_robotType == RobotType::MINI_CHEETAH) {
    balanceControllerVBL.set_inertia(0.025, 0.15, 0.18);
  }

  // Force limits
  threshold = this->_data->userParameters->two_leg_threshold;
  if (this->_data->_quadruped->_robotType == RobotType::CHEETAH_3)
    minForce = 10;
  else if (this->_data->_quadruped->_robotType == RobotType::MINI_CHEETAH)
    minForce = 1;
  maxForce = 70; // mass_in*9.81*1.6;
  for (int leg = 0; leg < 4; leg++) {
    minForces[leg] = minForce;
    maxForces[leg] = maxForce;
  }
}

/**
* Joint PD control to take desired legs out of contact
*/
template <typename T>
void FSM_State_TwoContactStand<T>::liftLeg(int leg, Vec3<T> qDes, Vec3<T> qdDes) {

  kpMat << 15, 0, 0, 0, 15, 0, 0, 0, 15;
  kdMat << 2, 0, 0, 0, 2, 0, 0, 0, 2;

  this->_data->_legController->commands[leg].kpJoint = kpMat;
  this->_data->_legController->commands[leg].kdJoint = kdMat;

  this->_data->_legController->commands[leg].qDes = qDes;
  this->_data->_legController->commands[leg].qdDes = qdDes;
}


/**
 * Gets the desired state of the robot from gamepad controls
 */
template <typename T>
void FSM_State_TwoContactStand<T>::get_desired_state() {
  // Nominal state
  pweight = 0.5;
  if (this->_data->_quadruped->_robotType == RobotType::CHEETAH_3)
    p_des[2] = 0.4;
  else if (this->_data->_quadruped->_robotType == RobotType::MINI_CHEETAH)
    p_des[2] = 0.225;
  for (int i = 0; i < 3; i++) {
    rpy[i] = 0.0;
    omegaDes[i] = 0.0;
    v_des[i] = 0.0;
  }

  // RC or gamepad to command desired state
  rpy[1] = 1.5 * this->_data->_desiredStateCommand->data.stateDes[4];
  rpy[2] = 0.35 * this->_data->_desiredStateCommand->data.stateDes[11];

  // Lift legs
  ramp_end = this->_data->userParameters->two_leg_ramp;
  if (this->_data->userParameters->stance_legs == 2) { // two legs
    pweight = 0.5 + 0.075 * this->_data->_desiredStateCommand->data.stateDes[6];
    p_des[0] = pweight * pFeet_yaw_world[3*1] + (1 - pweight) * pFeet_yaw_world[3*2];
    p_des[1] = pweight * pFeet_yaw_world[3*1+1] + (1 - pweight) * pFeet_yaw_world[3*2+1];

    if (contactStateScheduled[0] >= ramp_end)
      contactStateScheduled[0] -= 0.0002;
    else {
      contactStateScheduled[0] = 0.0;
      conPhase[0] = 0.0;
      if(pFeet_world_des[0 * 3 + 2] < 0.075)
        pFeet_world_des[0 * 3 + 2] += 0.001;
    }

    if (contactStateScheduled[3] >= ramp_end)
      contactStateScheduled[3] -= 0.0002;
    else {
      contactStateScheduled[3] = 0.0;
      conPhase[3] = 0.0;
      if(pFeet_world_des[3 * 3 + 2] < 0.075)
        pFeet_world_des[3 * 3 + 2] += 0.001;
    }

    q_lift_leg << 0., -1.45, 2.9;
  
  } else if (this->_data->userParameters->stance_legs == 3) { // three legs
    pweight = 0.425;
    p_des[0] = pweight * pFeet_yaw_world[3*0] + (1 - pweight) * pFeet_yaw_world[3*3];
    p_des[1] = pweight * pFeet_yaw_world[3*0+1] + (1 - pweight) * pFeet_yaw_world[3*3+1];

    if (contactStateScheduled[0] >= ramp_end)
      contactStateScheduled[0] -= 0.0002;
    else {
      contactStateScheduled[0] = 0.0;
      conPhase[0] = 0.0;
      if(pFeet_world_des[0 * 3 + 2] < 0.075)
        pFeet_world_des[0 * 3 + 2] += 0.001;
    }
      
    q_lift_leg << 0., -1.45, 2.9;

  } else { // four leg
    pweight = 0.5 + 0.075 * this->_data->_desiredStateCommand->data.stateDes[6];
    p_des[0] = pweight * pFeet_yaw_world[3*1] + (1 - pweight) * pFeet_yaw_world[3*2];
    pweight = 0.5 + 0.075 * this->_data->_desiredStateCommand->data.stateDes[7];
    p_des[1] = pweight * pFeet_yaw_world[3*1+1] + (1 - pweight) * pFeet_yaw_world[3*2+1];

    for (int i = 0; i < 4; i++){
      if(contactStateScheduled[i] < threshold){
        q_lift_leg << 0., -0.8, 2.15;
        double error;
        error = fabs(q_lift_leg[1]-this->_data->_legController->datas[i].q(1));
        if (error < 0.06) {
          contactStateScheduled[i] = 1;
          conPhase[i] = 0.5;
        }
      }
    }

  }

}

/**
 * Use data from leg controller and state estimator to compute positions
 * of the feet relative to the CoM
 */
template <typename T>
void FSM_State_TwoContactStand<T>::get_foot_locations() {

  for (int leg = 0; leg < 4; leg++) {
    
    // Compute vector from hip to foot (body coords)
    computeLegJacobianAndPosition(**&this->_data->_quadruped, this->_data->_legController->datas[leg].q,
                                  &this->_data->_legController->datas[leg].J, &pFeetVec, leg);
    
    // Following vectors express in world coords accounting for initial yaw offset
    pFeetVecBody = rBody_yaw * (this->_data->_quadruped->getHipLocation(leg) + pFeetVec); // Origin of body frame to foot

    pFeet_yaw[leg * 3] = (double)pFeetVecBody[0]-c_world[0]; // COM to foot
    pFeet_yaw[leg * 3 + 1] = (double)pFeetVecBody[1]-c_world[1];
    pFeet_yaw[leg * 3 + 2] = (double)pFeetVecBody[2]-c_world[2];

    pFeet_yaw_des[leg * 3] = (double)pFeetVecBody[0]-c_world[0]+p_act[0]-p_des[0]; // Desired COM to foot
    pFeet_yaw_des[leg * 3 + 1] = (double)pFeetVecBody[1]-c_world[1]+p_act[1]-p_des[1];
    pFeet_yaw_des[leg * 3 + 2] = (double)pFeetVecBody[2]-c_world[2]+p_act[2]-p_des[2];

    pFeet_yaw_world[leg * 3] = (double)p_act[0]+pFeet_yaw[leg * 3]; // absolute locations of the feet
    pFeet_yaw_world[leg * 3 + 1] = (double)p_act[1]+pFeet_yaw[leg * 3 + 1];
    pFeet_yaw_world[leg * 3 + 2] = (double)p_act[2]+pFeet_yaw[leg * 3 + 2];

    // Following vectors are expressed purely in the world frame (no yaw compensation - use for WBC)
    pFeetVecBody = this->_data->_stateEstimator->getResult().rBody.transpose() * (this->_data->_quadruped->getHipLocation(leg) + pFeetVec); // Origin of body frame to foot

    pFeet[leg * 3] = (double)pFeetVecBody[0]-c_world[0]; // COM to foot
    pFeet[leg * 3 + 1] = (double)pFeetVecBody[1]-c_world[1];
    pFeet[leg * 3 + 2] = (double)pFeetVecBody[2]-c_world[2];

    pFeet_world[leg * 3] = (double)p_act[0]+pFeet[leg * 3]; // absolute locations of the feet
    pFeet_world[leg * 3 + 1] = (double)p_act[1]+pFeet[leg * 3 + 1];
    pFeet_world[leg * 3 + 2] = (double)p_act[2]+pFeet[leg * 3 + 2];

    if(iter < 10){
      pFeet_world_des[leg * 3] = pFeet_world[leg * 3]; // Initial absolute position of the feet
      pFeet_world_des[leg * 3 + 1] = pFeet_world[leg * 3 + 1];
      pFeet_world_des[leg * 3 + 2] = pFeet_world[leg * 3 + 2];
    }

  }

}

/**
 * Visual debugging of coordinate frames and foot locations
 */
template <typename T>
void FSM_State_TwoContactStand<T>::debug_visualization() {

    // Draw the desired foot locations
    auto* FRd_Sphere = this->_data->visualizationData->addSphere();
    auto* FLd_Sphere = this->_data->visualizationData->addSphere();
    auto* BRd_Sphere = this->_data->visualizationData->addSphere();
    auto* BLd_Sphere = this->_data->visualizationData->addSphere();
    FRd_Sphere->position = pFoot_des[0];
    FRd_Sphere->radius = 0.02;
    FRd_Sphere->color = {0.6, 0.6, 0.2, 0.7};

    FLd_Sphere->position = pFoot_des[1];
    FLd_Sphere->radius = 0.02;
    FLd_Sphere->color = {0.6, 0.6, 0.2, 0.7};

    BRd_Sphere->position = pFoot_des[2];
    BRd_Sphere->radius = 0.02;
    BRd_Sphere->color = {0.6, 0.6, 0.2, 0.7};

    BLd_Sphere->position = pFoot_des[3];
    BLd_Sphere->radius = 0.02;
    BLd_Sphere->color = {0.6, 0.6, 0.2, 0.7};

    // Draw the world frame
    auto* worldX = this->_data->visualizationData->addArrow();
    auto* worldY = this->_data->visualizationData->addArrow();
    auto* worldZ = this->_data->visualizationData->addArrow();
    worldX->base_position = {0.0, 0.0, 0.0};
    worldX->direction = {0.5, 0.0, 0.0};
    worldX->color = {1.0, 0.0, 0.0, 0.7};
    worldX->shaft_width = 0.01;
    worldX->head_width = 0.02;
    worldX->head_length = 0.1;

    worldY->base_position = {0.0, 0.0, 0.0};
    worldY->direction = {0., 0.5, 0.0};
    worldY->color = {1.0, 0.0, 0.0, 0.7};
    worldY->shaft_width = 0.01;
    worldY->head_width = 0.02;
    worldY->head_length = 0.1;

    worldZ->base_position = {0.0, 0.0, 0.0};
    worldZ->direction = {0.0, 0.0, 0.5};
    worldZ->color = {1.0, 0.0, 0.0, 0.7};
    worldZ->shaft_width = 0.01;
    worldZ->head_width = 0.02;
    worldZ->head_length = 0.1;

    // Draw the body frame
    Vec3<float> unitX = {0.5,0.0,0.0};
    Vec3<float> unitY = {0.0,0.5,0.0};
    Vec3<float> unitZ = {0.0,0.0,0.5};
    auto* bodyX = this->_data->visualizationData->addArrow();
    auto* bodyY = this->_data->visualizationData->addArrow();
    auto* bodyZ = this->_data->visualizationData->addArrow();
    bodyX->base_position = {(float) p_act[0],(float) p_act[1],(float) p_act[2]};
    bodyX->direction = this->_data->_stateEstimator->getResult().rBody.transpose()*unitX;
    bodyX->color = {0.0, 0.0, 1.0, 0.7};
    bodyX->shaft_width = 0.01;
    bodyX->head_width = 0.02;
    bodyX->head_length = 0.1;

    bodyY->base_position = {(float)p_act[0],(float) p_act[1],(float) p_act[2]};
    bodyY->direction = this->_data->_stateEstimator->getResult().rBody.transpose()*unitY;
    bodyY->color = {0.0, 0.0, 1.0, 0.7};
    bodyY->shaft_width = 0.01;
    bodyY->head_width = 0.02;
    bodyY->head_length = 0.1;

    bodyZ->base_position = {(float)p_act[0],(float) p_act[1],(float) p_act[2]};
    bodyZ->direction = this->_data->_stateEstimator->getResult().rBody.transpose()*unitZ;
    bodyZ->color = {0.0, 0.0, 1.0, 0.7};
    bodyZ->shaft_width = 0.01;
    bodyZ->head_width = 0.02;
    bodyZ->head_length = 0.1;

}

 /**
 * Convert euler angles to rotation matrix
 */
template <typename T>
void FSM_State_TwoContactStand<T>::rpyToR(Mat3<float> &R, double* rpy_in) {

   Mat3<float> Rz, Ry, Rx;

   Rz.setIdentity();
   Ry.setIdentity();
   Rx.setIdentity();

   Rz(0,0) = cos(rpy_in[2]);
   Rz(0,1) = -sin(rpy_in[2]);
   Rz(1,0) = sin(rpy_in[2]);
   Rz(1,1) = cos(rpy_in[2]);

   Ry(0,0) = cos(rpy_in[1]);
   Ry(0,2) = sin(rpy_in[1]);
   Ry(2,0) = -sin(rpy_in[1]);
   Ry(2,2) = cos(rpy_in[1]);

   Rx(1,1) = cos(rpy_in[0]);
   Rx(1,2) = -sin(rpy_in[0]);
   Rx(2,1) = sin(rpy_in[0]);
   Rx(2,2) = cos(rpy_in[0]);

   R = Rz*Ry*Rx;

 }

/**
 * Convert euler angles to quaternions
 */
template <typename T>
void FSM_State_TwoContactStand<T>::eulerToQuat(double* rpy_in, double* quat_in) {

  quat_in[0] = cos(rpy_in[2]*0.5)*cos(rpy_in[1]*0.5)*cos(rpy_in[0]*0.5) + sin(rpy_in[2]*0.5)*sin(rpy_in[1]*0.5)*sin(rpy_in[0]*0.5);
  quat_in[1] = cos(rpy_in[2]*0.5)*cos(rpy_in[1]*0.5)*sin(rpy_in[0]*0.5) - sin(rpy_in[2]*0.5)*sin(rpy_in[1]*0.5)*cos(rpy_in[0]*0.5);
  quat_in[2] = sin(rpy_in[2]*0.5)*cos(rpy_in[1]*0.5)*sin(rpy_in[0]*0.5) + cos(rpy_in[2]*0.5)*sin(rpy_in[1]*0.5)*cos(rpy_in[0]*0.5);
  quat_in[3] = sin(rpy_in[2]*0.5)*cos(rpy_in[1]*0.5)*cos(rpy_in[0]*0.5) - cos(rpy_in[2]*0.5)*sin(rpy_in[1]*0.5)*sin(rpy_in[0]*0.5);

}

/**
 * Convert quaternions to euler angles
 */
template <typename T>
void FSM_State_TwoContactStand<T>::quatToEuler(double* quat_in, double* rpy_in) {
  // roll (x-axis rotation)
  double sinr_cosp = +2.0 * (quat_in[0] * quat_in[1] + quat_in[2] * quat_in[3]);
  double cosr_cosp = +1.0 - 2.0 * (quat_in[1] * quat_in[1] + quat_in[2] * quat_in[2]);
  rpy_in[0] = atan2(sinr_cosp, cosr_cosp);

  // pitch (y-axis rotation)
  double sinp = +2.0 * (quat_in[0] * quat_in[2] - quat_in[3] * quat_in[1]);
  if (fabs(sinp) >= 1)
    rpy_in[1] = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
  else
    rpy_in[1] = asin(sinp);

  // yaw (z-axis rotation)
  double siny_cosp = +2.0 * (quat_in[0] * quat_in[3] + quat_in[1] * quat_in[2]);
  double cosy_cosp = +1.0 - 2.0 * (quat_in[2] * quat_in[2] + quat_in[3] * quat_in[3]);  
  rpy_in[2] = atan2(siny_cosp, cosy_cosp);
}


/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_TwoContactStand<T>::checkTransition() {
  this->nextStateName = this->stateName;
  iter++;

  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
    case K_TWO_CONTACT_STAND:
      break;

    case K_RECOVERY_STAND:
      this->nextStateName = FSM_StateName::RECOVERY_STAND;
      break;

    case K_LOCOMOTION:
      this->nextStateName = FSM_StateName::LOCOMOTION;
      break;

    case K_BOUNDING:
      this->nextStateName = FSM_StateName::BOUNDING;
      break;


    case K_PASSIVE:  // normal c
      this->nextStateName = FSM_StateName::PASSIVE;
      break;

    case K_BALANCE_STAND: 
      this->nextStateName = FSM_StateName::BALANCE_STAND;
      break;

    default:
      std::cout << "[CONTROL FSM] Bad Request: Cannot transition from "
                << K_TWO_CONTACT_STAND << " to "
                << this->_data->controlParameters->control_mode << std::endl;
  }

  // Get the next state
  return this->nextStateName;
}

/**
 * Handles the actual transition for the robot between states.
 * Returns true when the transition is completed.
 *
 * @return true if transition is complete
 */
template <typename T>
TransitionData<T> FSM_State_TwoContactStand<T>::transition() {
  // Finish Transition
  switch (this->nextStateName) {
    case FSM_StateName::PASSIVE:  // normal
      this->transitionData.done = true;
      break;

    case FSM_StateName::BALANCE_STAND:
      this->transitionData.done = true;
      break;

    case FSM_StateName::LOCOMOTION:
      this->transitionData.done = true;
      break;

    case FSM_StateName::BOUNDING:
      this->transitionData.done = true;
      break;

    case FSM_StateName::RECOVERY_STAND:
      this->transitionData.done = true;
      break;

    default:
      std::cout << "[CONTROL FSM] Something went wrong in transition"
                << std::endl;
  }

  // Return the transition data to the FSM
  return this->transitionData;
}

/**
 * Cleans up the state information on exiting the state.
 */
template <typename T>
void FSM_State_TwoContactStand<T>::onExit() {
  // Nothing to clean up when exiting
}



// template class FSM_State_TwoContactStand<double>;
template class FSM_State_TwoContactStand<float>;
