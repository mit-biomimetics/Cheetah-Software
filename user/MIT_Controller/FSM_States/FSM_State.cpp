/*============================= FSM State =============================*/
/**
 * FSM State base class
 */

#include "FSM_State.h"

/**
 * Constructor for the FSM State class.
 *
 * @param _controlFSMData holds all of the relevant control data
 * @param stateNameIn the enumerated state name
 * @param stateStringIn the string name of the current FSM state
 */
template <typename T>
FSM_State<T>::FSM_State(ControlFSMData<T>* _controlFSMData,
                        FSM_StateName stateNameIn, std::string stateStringIn)
    : _data(_controlFSMData),
      stateName(stateNameIn),
      stateString(stateStringIn) {
  transitionData.zero();
  std::cout << "[FSM_State] Initialized FSM state: " << stateStringIn
            << std::endl;
}

/**
 * Cartesian impedance control for a given leg.
 *
 * @param leg the leg number to control
 * @param qDes desired joint position
 * @param dqDes desired joint velocity
 */
template <typename T>
void FSM_State<T>::jointPDControl(
    int leg, Vec3<T> qDes, Vec3<T> qdDes) {

  kpMat << 80, 0, 0, 0, 80, 0, 0, 0, 80;
  kdMat << 1, 0, 0, 0, 1, 0, 0, 0, 1;

  _data->_legController->commands[leg].kpJoint = kpMat;
  _data->_legController->commands[leg].kdJoint = kdMat;

  _data->_legController->commands[leg].qDes = qDes;
  _data->_legController->commands[leg].qdDes = qdDes;
}

/**
 * Cartesian impedance control for a given leg.
 *
 * @param leg the leg number to control
 * @param pDes desired foot position
 * @param vDes desired foot velocity
 * @param kp_cartesian P gains
 * @param kd_cartesian D gains
 */
template <typename T>
void FSM_State<T>::cartesianImpedanceControl(int leg, Vec3<T> pDes,
                                             Vec3<T> vDes,
                                             Vec3<double> kp_cartesian,
                                             Vec3<double> kd_cartesian) {
  _data->_legController->commands[leg].pDes = pDes;
  // Create the cartesian P gain matrix
  kpMat << kp_cartesian[0], 0, 0, 0,
      kp_cartesian[1], 0, 0, 0,
      kp_cartesian[2];
  _data->_legController->commands[leg].kpCartesian = kpMat;

  _data->_legController->commands[leg].vDes = vDes;
  // Create the cartesian D gain matrix
  kdMat << kd_cartesian[0], 0, 0, 0, kd_cartesian[1], 0, 0, 0, kd_cartesian[2];
  _data->_legController->commands[leg].kdCartesian = kdMat;
}

/**
 *
 */
/*
template <typename T>
void FSM_State<T>::footstepHeuristicPlacement(int leg) {

  // Create the projection matrix for the 2D foot placement components
  Mat23<float> projectionMatrix;
  projectionMatrix << 1, 0, 0,
                   0, 1, 0;

  Vec3<float> velDes = _data->_desiredStateCommand->data.stateDes.block<3, 1>(6,
0); Vec3<float> angVelDes = _data->_desiredStateCommand->data.stateDes.block<3,
1>(9, 0); Mat3<float> rBody = _data->_stateEstimate.rBody;

  // Find each of the footstep locations for the swing feet
  for (int leg = 0; leg < 4; leg++) {
    if (_data->_gaitScheduler->gaitData.contactStateScheduled(leg)) {
      // The leg is in contact so nothing to do here

    } else {
      if (_data->_gaitScheduler->gaitData._currentGait ==
GaitType::TRANSITION_TO_STAND) {
        // Position the legs under the hips to stand...
        // Could also get rid of this and simply send 0 velocity ang vel
commands
        // from the CoM desired planner...
        Vec3<float> posHip = _data->_quadruped.getHipLocation(leg);
        footstepLocations.col(leg) <<
projectionMatrix.transpose()*projectionMatrix*
                                   (_data->_stateEstimate.position + // rBody *
posHip); } else {
        // Pull out the approximate yaw rate component of the robot in the
world. Vec3<float> yaw_rate; yaw_rate << 0, 0, _stateEstimate.omegaWorld(3);

        Vec3<float> posHip = _data->_quadruped.getHipLocation(leg);

        float timeStance = _data->_gaitScheduler->gaitData.timeStance(leg);

        // Footstep heuristic composed of several parts in the world frame
        footstepLocations.col(leg) <<
projectionMatrix.transpose()*projectionMatrix*      // Ground projection
                                   (_stateEstimate.position + // rBody * posHip
+                                      // Foot under hips timeStance / 2 *
velDes +                             // Raibert Heuristic timeStance / 2 *
(angVelDes.cross(rBody * posHip)) +  // Turning Raibert Heuristic
                                    (_stateEstimate.vBody - velDes));
      }
    }

  }
}
*/

/**
 * Gait independent formulation for choosing appropriate GRF and step locations
 * as well as converting them to leg controller understandable values.
 */
template <typename T>
void FSM_State<T>::runControls() {
  // This option should be set from the user interface or autonomously
  // eventually
  int CONTROLLER_OPTION = 1;

  // Reset the forces and steps to 0
  footFeedForwardForces = Mat34<T>::Zero();
  footstepLocations = Mat34<T>::Zero();

  // Choose the controller to run for picking step locations and balance forces
  if (CONTROLLER_OPTION == 0) {
    // Test to make sure we can control the robot
    // Test to make sure we can control the robot these will be calculated by
    // the controllers
    for (int leg = 0; leg < 4; leg++) {
      footFeedForwardForces.col(leg) << 0.0, 0.0, 0;  //-220.36;
      // footFeedForwardForces.col(leg) = stateEstimate.rBody *
      // footFeedForwardForces.col(leg);

      footstepLocations.col(leg) << 0.0, 0.0,
          -_data->_quadruped->_maxLegLength / 2;
    }
  } else if (CONTROLLER_OPTION == 1) {
    // QP Balance Controller
    runBalanceController();

    // Swing foot landing positions are calculated with heuristics
    for (int leg = 0; leg < 4; leg++) {
      footstepLocations.col(leg) << 0.0, 0.0,
          -_data->_quadruped->_maxLegLength / 2;
    }  // footstepHeuristicPlacement();

  } else if (CONTROLLER_OPTION == 2) {
    // WBC
    runWholeBodyController();

  } else if (CONTROLLER_OPTION == 3) {
    // cMPC
    runConvexModelPredictiveController();

    // Swing foot landing positions are calculated with heuristics
    // footstepHeuristicPlacement();

  } else if (CONTROLLER_OPTION == 4) {
    // RPC
    runRegularizedPredictiveController();

  } else {
    // Zero out the commands if a controller was not selected
    jointFeedForwardTorques =
        Mat34<float>::Zero();                // feed forward joint torques
    jointPositions = Mat34<float>::Zero();   // joint angle positions
    jointVelocities = Mat34<float>::Zero();  // joint angular velocities
    footFeedForwardForces =
        Mat34<float>::Zero();              // feedforward forces at the feet
    footPositions = Mat34<float>::Zero();  // cartesian foot positions
    footVelocities = Mat34<float>::Zero();

    // Print an error message
    std::cout << "[FSM_State] ERROR: No known controller was selected: "
              << CONTROLLER_OPTION << std::endl;
  }
}

/**
 *
 */
template <typename T>
void FSM_State<T>::runBalanceController() {
  double minForce = 25;
  double maxForce = 500;
  double contactStateScheduled[4];  // = {1, 1, 1, 1};
  for (int i = 0; i < 4; i++) {
    contactStateScheduled[i] =
        _data->_gaitScheduler->gaitData.contactStateScheduled(i);
  }

  double minForces[4];  // = {minForce, minForce, minForce, minForce};
  double maxForces[4];  // = {maxForce, maxForce, maxForce, maxForce};
  for (int leg = 0; leg < 4; leg++) {
    minForces[leg] = contactStateScheduled[leg] * minForce;
    maxForces[leg] = contactStateScheduled[leg] * maxForce;
  }

  double COM_weights_stance[3] = {1, 1, 10};
  double Base_weights_stance[3] = {20, 10, 10};
  double pFeet[12], p_des[3], p_act[3], v_des[3], v_act[3], O_err[3], rpy[3],
      omegaDes[3];
  double se_xfb[13];
  double kpCOM[3], kdCOM[3], kpBase[3], kdBase[3];

  for (int i = 0; i < 4; i++) {
    se_xfb[i] = (double)_data->_stateEstimator->getResult().orientation(i);
  }
  // se_xfb[3] = 1.0;
  for (int i = 0; i < 3; i++) {
    rpy[i] = 0;  //(double)_data->_stateEstimator->getResult().rpy(i);
    p_des[i] = (double)_data->_stateEstimator->getResult().position(i);
    p_act[i] = (double)_data->_stateEstimator->getResult().position(i);
    omegaDes[i] =
        0;  //(double)_data->_stateEstimator->getResult().omegaBody(i);
    v_act[i] = (double)_data->_stateEstimator->getResult().vBody(i);
    v_des[i] = (double)_data->_stateEstimator->getResult().vBody(i);

    se_xfb[4 + i] = (double)_data->_stateEstimator->getResult().position(i);
    se_xfb[7 + i] = (double)_data->_stateEstimator->getResult().omegaBody(i);
    se_xfb[10 + i] = (double)_data->_stateEstimator->getResult().vBody(i);

    // Set the translational and orientation gains
    kpCOM[i] = (double)_data->controlParameters->kpCOM(i);
    kdCOM[i] = (double)_data->controlParameters->kdCOM(i);
    kpBase[i] = (double)_data->controlParameters->kpBase(i);
    kdBase[i] = (double)_data->controlParameters->kdBase(i);
  }

  Vec3<T> pFeetVec;
  Vec3<T> pFeetVecCOM;
  // Get the foot locations relative to COM
  for (int leg = 0; leg < 4; leg++) {
    computeLegJacobianAndPosition(**&_data->_quadruped,
                                  _data->_legController->datas[leg].q,
                                  (Mat3<T>*)nullptr, &pFeetVec, 1);
    //pFeetVecCOM = _data->_stateEstimator->getResult().rBody.transpose() *
                  //(_data->_quadruped->getHipLocation(leg) + pFeetVec);

    pFeetVecCOM = _data->_stateEstimator->getResult().rBody.transpose() *
                  (_data->_quadruped->getHipLocation(leg) + _data->_legController->datas[leg].p);


    pFeet[leg * 3] = (double)pFeetVecCOM[0];
    pFeet[leg * 3 + 1] = (double)pFeetVecCOM[1];
    pFeet[leg * 3 + 2] = (double)pFeetVecCOM[2];
  }
  
  balanceController.set_alpha_control(0.01);
  balanceController.set_friction(0.5);
  balanceController.set_mass(46.0);
  balanceController.set_wrench_weights(COM_weights_stance, Base_weights_stance);
  balanceController.set_PDgains(kpCOM, kdCOM, kpBase, kdBase);
  balanceController.set_desiredTrajectoryData(rpy, p_des, omegaDes, v_des);
  balanceController.SetContactData(contactStateScheduled, minForces, maxForces);
  balanceController.updateProblemData(se_xfb, pFeet, p_des, p_act, v_des, v_act,
                                      O_err, 0.0);

  double fOpt[12];
  balanceController.solveQP_nonThreaded(fOpt);

  // Publish the results over LCM
  balanceController.publish_data_lcm();

  // Copy the results to the feed forward forces
  for (int leg = 0; leg < 4; leg++) {
    footFeedForwardForces.col(leg) << (T)fOpt[leg * 3], (T)fOpt[leg * 3 + 1],
        (T)fOpt[leg * 3 + 2];
  }
}

/**
 * Gait independent formulation for choosing appropriate GRF and step locations
 * as well as converting them to leg controller understandable values.
 */
template <typename T>
void FSM_State<T>::turnOnAllSafetyChecks() {
  // Pre controls safety checks
  checkSafeOrientation = true;  // check roll and pitch

  // Post control safety checks
  checkPDesFoot = true;          // do not command footsetps too far
  checkForceFeedForward = true;  // do not command huge forces
  checkLegSingularity = true;    // do not let leg
}

/**
 *
 */
template <typename T>
void FSM_State<T>::turnOffAllSafetyChecks() {
  // Pre controls safety checks
  checkSafeOrientation = false;  // check roll and pitch

  // Post control safety checks
  checkPDesFoot = false;          // do not command footsetps too far
  checkForceFeedForward = false;  // do not command huge forces
  checkLegSingularity = false;    // do not let leg
}

// template class FSM_State<double>;
template class FSM_State<float>;
