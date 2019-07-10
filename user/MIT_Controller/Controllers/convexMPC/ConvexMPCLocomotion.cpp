#include <iostream>
#include <Utilities/Timer.h>

#include "ConvexMPCLocomotion.h"
#include "convexMPC_interface.h"

#define DRAW_DEBUG_SWINGS
//#define DRAW_DEBUG_PATH


///////////////
// GAIT
///////////////
Gait::Gait(int nMPC_segments, Vec4<int> offsets, Vec4<int> durations, const std::string &name) :
  _offsets(offsets.array()),
  _durations(durations.array()),
  _nIterations(nMPC_segments)
{
  _mpc_table = new int[nMPC_segments * 4];

  _offsetsFloat = offsets.cast<float>() / (float) nMPC_segments;
  _durationsFloat = durations.cast<float>() / (float) nMPC_segments;
  std::cout << "Gait " << name << "\n";
  std::cout << "nMPC_segments    : " << _nIterations << "\n";
  std::cout << "offsets (int)    : " << _offsets.transpose() << "\n";
  std::cout << "durations (int)  : " << _durations.transpose() << "\n";
  std::cout << "offsets (float)  : " << _offsetsFloat.transpose() << "\n";
  std::cout << "durations (float): " << _durationsFloat.transpose() << "\n";
  std::cout << "\n\n";

  _stance = durations[0];
  _swing = nMPC_segments - durations[0];

}


Gait::~Gait() {
  delete[] _mpc_table;
}


Vec4<float> Gait::getContactState()
{
  Array4f progress = _phase - _offsetsFloat;

  for(int i = 0; i < 4; i++)
  {
    if(progress[i] < 0) progress[i] += 1.;
    if(progress[i] > _durationsFloat[i])
    {
      progress[i] = 0.;
    }
    else
    {
      progress[i] = progress[i] / _durationsFloat[i];
    }
  }

  return progress.matrix();
}

Vec4<float> Gait::getSwingState()
{
  Array4f swing_offset = _offsetsFloat + _durationsFloat;
  for(int i = 0; i < 4; i++)
    if(swing_offset[i] > 1) swing_offset[i] -= 1.;
  Array4f swing_duration = 1. - _durationsFloat;

  Array4f progress = _phase - swing_offset;

  for(int i = 0; i < 4; i++)
  {
    if(progress[i] < 0) progress[i] += 1.f;
    if(progress[i] > swing_duration[i])
    {
      progress[i] = 0.;
    }
    else
    {
      progress[i] = progress[i] / swing_duration[i];
    }
  }

  return progress.matrix();
}

int* Gait::mpc_gait()
{
  for(int i = 0; i < _nIterations; i++)
  {
    int iter = (i + _iteration + 1) % _nIterations;
    Array4i progress = iter - _offsets;
    for(int j = 0; j < 4; j++)
    {
      if(progress[j] < 0) progress[j] += _nIterations;
      if(progress[j] < _durations[j])
        _mpc_table[i*4 + j] = 1;
      else
        _mpc_table[i*4 + j] = 0;
    }
  }

  return _mpc_table;
}

void Gait::setIterations(int iterationsPerMPC, int currentIteration)
{
  _iteration = (currentIteration / iterationsPerMPC) % _nIterations;
  _phase = (float)(currentIteration % (iterationsPerMPC * _nIterations)) / (float) (iterationsPerMPC * _nIterations);
}


////////////////////
// Controller
////////////////////

ConvexMPCLocomotion::ConvexMPCLocomotion() :
  horizonLength(10),
  trotting(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(5,5,5,5),"Trotting"),
  bounding(horizonLength, Vec4<int>(5,5,0,0),Vec4<int>(5,5,5,5),"Bounding"),
  pronking(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(4,4,4,4),"Pronking"),
  galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(6,6,6,6),"Galloping"),
  standing(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(10,10,10,10),"Standing"),
  trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(3,3,3,3),"Trot Running"),
  walking(horizonLength, Vec4<int>(0,3,5,8), Vec4<int>(5,5,5,5), "Walking"),
  walking2(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(7,7,7,7), "Walking2"),
  pacing(horizonLength, Vec4<int>(5,0,5,0),Vec4<int>(5,5,5,5),"Pacing")
{
  dtMPC = 0.001 * iterationsBetweenMPC;
  setup_problem(dtMPC, horizonLength, 0.4, 120);
  rpy_comp[0] = 0;
  rpy_comp[1] = 0;
  rpy_comp[2] = 0;
  rpy_int[0] = 0;
  rpy_int[1] = 0;
  rpy_int[2] = 0;

  for(int i = 0; i < 4; i++)
    firstSwing[i] = true;
}

template<>
void ConvexMPCLocomotion::run(ControlFSMData<float>& data) {
  bool omniMode = false;
  gaitNumber = data.userParameters->cmpc_gait;
  if(gaitNumber >= 10) {
    gaitNumber -= 10;
    omniMode = true;
  }

//  auto* debugSphere = data.visualizationData->addSphere();
//  debugSphere->color = {1,1,1,0.5};
//  debugSphere->radius = 1;

  auto& seResult = data._stateEstimator->getResult();
  auto& stateCommand = data._desiredStateCommand;

  // Check if transition to standing
  if(((gaitNumber == 4) && current_gait != 4) || firstRun)
  {
    printf("Transition to standing\n");
    stand_traj[0] = seResult.position[0];
    stand_traj[1] = seResult.position[1];
    stand_traj[2] = 0.21;
    stand_traj[3] = 0;
    stand_traj[4] = 0;
    stand_traj[5] = seResult.rpy[2];
    world_position_desired[0] = stand_traj[0];
    world_position_desired[1] = stand_traj[1];
  }



  // pick gait
  Gait* gait = &trotting;
  if(gaitNumber == 1)
    gait = &bounding;
  else if(gaitNumber == 2)
    gait = &pronking;
  else if(gaitNumber == 3)
    gait = &galloping;
  else if(gaitNumber == 4)
    gait = &standing;
  else if(gaitNumber == 5)
    gait = &trotRunning;
  else if(gaitNumber == 6)
    gait = &walking;
  else if(gaitNumber == 7)
    gait = &walking2;
  else if(gaitNumber == 8)
    gait = &pacing;
  current_gait = gaitNumber;

  // integrate position setpoint
  Vec3<float> v_des_robot(stateCommand->data.stateDes[6], stateCommand->data.stateDes[7],0);
  Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
  Vec3<float> v_robot = seResult.vWorld;


  //Integral-esque pitche and roll compensation
  if(fabs(v_robot[0]) > .2)   //avoid dividing by zero
  {
    rpy_int[1] += 0.001*(stateCommand->data.stateDes[4] /*-hw_i->state_estimator->se_ground_pitch*/ - seResult.rpy[1])/v_robot[0];
  }
  if(fabs(v_robot[1]) > 0.1)
  {
    rpy_int[0] += 0.001*(stateCommand->data.stateDes[3] /*-hw_i->state_estimator->se_ground_pitch*/ - seResult.rpy[0])/v_robot[1];
  }

  rpy_int[0] = fminf(fmaxf(rpy_int[0], -.25), .25);
  rpy_int[1] = fminf(fmaxf(rpy_int[1], -.25), .25);
  rpy_comp[1] = v_robot[0] * rpy_int[1];
  rpy_comp[0] = v_robot[1] * rpy_int[0] * (gaitNumber!=8);  //turn off for pronking




  for(int i = 0; i < 4; i++) {
    pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
  }

  if(gait != &standing) {
    world_position_desired += .001 * Vec3<float>(v_des_world[0], v_des_world[1], 0);
  }

  // some first time initialization
  if(firstRun)
  {
    world_position_desired[0] = seResult.position[0];
    world_position_desired[1] = seResult.position[1];
    world_position_desired[2] = seResult.rpy[2];

    for(int i = 0; i < 4; i++)
    {

      footSwingTrajectories[i].setHeight(0.1);
      footSwingTrajectories[i].setInitialPosition(pFoot[i]);
      footSwingTrajectories[i].setFinalPosition(pFoot[i]);

    }
    firstRun = false;
  }

  // foot placement
  swingTimes[0] = dtMPC * gait->_swing;
  swingTimes[1] = dtMPC * gait->_swing;
  swingTimes[2] = dtMPC * gait->_swing;
  swingTimes[3] = dtMPC * gait->_swing;
  float side_sign[4] = {-1, 1, -1, 1};
  for(int i = 0; i < 4; i++)
  {

    if(firstSwing[i]) {
      swingTimeRemaining[i] = swingTimes[i];
    } else {
      swingTimeRemaining[i] -= 0.001f;
    }
    //if(firstSwing[i]) {
      footSwingTrajectories[i].setHeight(.1);
      Vec3<float> offset(0, side_sign[i] * .065, 0);

      Vec3<float> pRobotFrame = (data._quadruped->getHipLocation(i) + offset);
      Vec3<float> pYawCorrected = coordinateRotation(CoordinateAxis::Z, -stateCommand->data.stateDes[11] * gait->_stance * dtMPC / 2) * pRobotFrame;


      Vec3<float> Pf = seResult.position +
                       seResult.rBody.transpose() * pYawCorrected + seResult.vWorld * swingTimeRemaining[i];

      float p_rel_max = 0.3f;
      float pfx_rel = seResult.vWorld[0] * .5 * gait->_stance * dtMPC +
                      .03f*(seResult.vWorld[0]-v_des_world[0]) +
                      (0.5f*seResult.position[2]/9.81f) * (seResult.vWorld[1]*stateCommand->data.stateDes[2]);
      float pfy_rel = seResult.vWorld[1] * .5 * gait->_stance * dtMPC +
                      .03f*(seResult.vWorld[1]-v_des_world[1]) +
                      (0.5f*seResult.position[2]/9.81f) * (-seResult.vWorld[0]*stateCommand->data.stateDes[2]);
      pfx_rel = fminf(fmaxf(pfx_rel, -p_rel_max), p_rel_max);
      pfy_rel = fminf(fmaxf(pfy_rel, -p_rel_max), p_rel_max);
      Pf[0] +=  pfx_rel;
      Pf[1] +=  pfy_rel;
      Pf[2] = -0.01;
      footSwingTrajectories[i].setFinalPosition(Pf);
    //}

  }

  // swing timess

//  for(int i = 0; i < 4; i++)
//  {
//    footSwingController[i]->setSwingTime(swingTimes[i]);
//  }


  // calc gait
  gait->setIterations(iterationsBetweenMPC, iterationCounter);
  iterationCounter++;

  // load LCM leg swing gains
  Kp << 700, 0, 0,
    0, 700, 0,
    0, 0, 150;
  Kp_stance = 0*Kp;


  Kd << 11, 0, 0,
    0, 11, 0,
    0, 0, 11;
  Kd_stance = Kd;
  // gait
  Vec4<float> contactStates = gait->getContactState();
  Vec4<float> swingStates = gait->getSwingState();
  int* mpcTable = gait->mpc_gait();
  updateMPCIfNeeded(mpcTable, data, omniMode);

//  StateEstimator* se = hw_i->state_estimator;
  Vec4<float> se_contactState(0,0,0,0);

#ifdef DRAW_DEBUG_PATH
  auto* trajectoryDebug = data.visualizationData->addPath();
  if(trajectoryDebug) {
    trajectoryDebug->num_points = 10;
    trajectoryDebug->color = {0.2, 0.2, 0.7, 0.5};
    for(int i = 0; i < 10; i++) {
      trajectoryDebug->position[i][0] = trajAll[12*i + 3];
      trajectoryDebug->position[i][1] = trajAll[12*i + 4];
      trajectoryDebug->position[i][2] = trajAll[12*i + 5];
      auto* ball = data.visualizationData->addSphere();
      ball->radius = 0.01;
      ball->position = trajectoryDebug->position[i];
      ball->color = {1.0, 0.2, 0.2, 0.5};
    }
  }
#endif

  for(int foot = 0; foot < 4; foot++)
  {
    float contactState = contactStates[foot];
    float swingState = swingStates[foot];
    if(swingState > 0) // foot is in swing
    {
      if(firstSwing[foot])
      {
        firstSwing[foot] = false;
        footSwingTrajectories[foot].setInitialPosition(pFoot[foot]);
      }

#ifdef DRAW_DEBUG_SWINGS
      auto* debugPath = data.visualizationData->addPath();
      if(debugPath) {
        debugPath->num_points = 100;
        debugPath->color = {0.2,1,0.2,0.5};
        float step = (1.f - swingState) / 100.f;
        for(int i = 0; i < 100; i++) {
          footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState + i * step, swingTimes[foot]);
          debugPath->position[i] = footSwingTrajectories[foot].getPosition();
        }
      }
      auto* finalSphere = data.visualizationData->addSphere();
      if(finalSphere) {
        finalSphere->position = footSwingTrajectories[foot].getPosition();
        finalSphere->radius = 0.02;
        finalSphere->color = {0.6, 0.6, 0.2, 0.7};
      }
      footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState, swingTimes[foot]);
      auto* actualSphere = data.visualizationData->addSphere();
      auto* goalSphere = data.visualizationData->addSphere();
      goalSphere->position = footSwingTrajectories[foot].getPosition();
      actualSphere->position = pFoot[foot];
      goalSphere->radius = 0.02;
      actualSphere->radius = 0.02;
      goalSphere->color = {0.2, 1, 0.2, 0.7};
      actualSphere->color = {0.8, 0.2, 0.2, 0.7};
#endif
      footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState, swingTimes[foot]);


//      footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
//                                          hw_i->leg_controller->leg_datas[foot].qd, 0); // velocity dependent friction compensation todo removed
      //hw_i->leg_controller->leg_datas[foot].qd, fsm->main_control_settings.variable[2]);

      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);



//      Vec3<float> pFootBody = data._legController->datas[foot].p + data._quadruped->getHipLocation(foot);
//      for(int i = 0; i<4; i++) todo removed
//      {
//        if(i != foot)
//        {
//          vec3 pOtherFoot = hw_i->leg_controller->leg_datas[i].p + hw_i->state_estimator->hip_positions.col(i);
//          vec3 distanceVec = pOtherFoot-pFootBody;
//          distanceVec(2) = 0; //don't care about z
//          float distance = distanceVec.norm();
//          if(distance<0.05)
//          {
//            hw_i->leg_controller->leg_commands[foot].f_ff = -1000*(.05 - distance)*(distanceVec/distance);
//            cout << "Leg " << foot << " anti-swing collision with leg: " << i << " distance: " << distance << "\n";
//          }
//        }
//      }

      //cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";
      data._legController->commands[foot].pDes = pDesLeg;
      data._legController->commands[foot].vDes = vDesLeg;
      data._legController->commands[foot].kpCartesian = Kp;
      data._legController->commands[foot].kdCartesian = Kd;

      //singularity barrier
      data._legController->commands[foot].tauFeedForward[2] = 50*(data._legController->datas[foot].q(2)<.1)*data._legController->datas[foot].q(2);

      //Vec3<float> pErr = data._legController->datas[foot].p - data._legController->commands[foot].pDes;
      //printf("[%d] %.3f %.3f %.3f\n", foot, pErr[0], pErr[1], pErr[2]);
//            vec3 perr = pDesLeg - hw_i->leg_controller->leg_datas[foot].p;
//            vec3 verr = vDesLeg - hw_i->leg_controller->leg_datas[foot].v;
//            cout << "FOOT: " << foot << "\nperr: " << perr(2) << "\npdes: " << pDesLeg(2) << "\n";
      //cout << "Tau ff: " << footSwingController[foot]->getTauFF().transpose() << "\n";
    }
    else // foot is in stance
    {
      firstSwing[foot] = true;

#ifdef DRAW_DEBUG_SWINGS
      auto* actualSphere = data.visualizationData->addSphere();
      actualSphere->position = pFoot[foot];
      actualSphere->radius = 0.02;
      actualSphere->color = {0.2, 0.2, 0.8, 0.7};
#endif

      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);
      //cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";
      data._legController->commands[foot].pDes = pDesLeg;
      data._legController->commands[foot].vDes = vDesLeg;
      data._legController->commands[foot].kpCartesian = Kp_stance;
      data._legController->commands[foot].kdCartesian = Kd_stance;

      data._legController->commands[foot].forceFeedForward = f_ff[foot];
      data._legController->commands[foot].kdJoint = Mat3<float>::Identity() * 0.2;

//      footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
//                                          hw_i->leg_controller->leg_datas[foot].qd, 0); todo removed
      // hw_i->leg_controller->leg_commands[foot].tau_ff += 0*footSwingController[foot]->getTauFF();
      data._legController->commands[foot].tauFeedForward[2] = 50*(data._legController->datas[foot].q(2)<.1)*data._legController->datas[foot].q(2);
//            cout << "Foot " << foot << " force: " << f_ff[foot].transpose() << "\n";
      se_contactState[foot] = contactState;
    }
  }

  // se->set_contact_state(se_contactState); todo removed
  data._stateEstimator->setContactPhase(se_contactState);

}

template<>
void ConvexMPCLocomotion::run(ControlFSMData<double>& data) {
  (void)data;
  printf("call to old CMPC with double!\n");

}

void ConvexMPCLocomotion::updateMPCIfNeeded(int *mpcTable, ControlFSMData<float> &data, bool omniMode) {
  iterationsBetweenMPC = 30;
  if((iterationCounter % iterationsBetweenMPC) == 0)
  {
    auto seResult = data._stateEstimator->getResult();
    auto& stateCommand = data._desiredStateCommand;
    float* p = seResult.position.data();
    float* v = seResult.vWorld.data();
    float* w = seResult.omegaWorld.data();
    float* q = seResult.orientation.data();
    float r[12];
    for(int i = 0; i < 12; i++)
      r[i] = pFoot[i%4][i/4]  - seResult.position[i/4];

    float Q[12] = {0.25, 0.25, 10, 2, 2, 20, 0, 0, 0.3, 0.2, 0.2, 0.2};
    float yaw = seResult.rpy[2];
    float* weights = Q;
    float alpha = 4e-5; // make setting eventually

    //printf("current posistion: %3.f %.3f %.3f\n", p[0], p[1], p[2]);

    if(alpha > 1e-4)
    {

      std::cout << "Alpha was set too high (" << alpha << ") adjust to 1e-5\n";
      alpha = 1e-5;
    }
    Vec3<float> v_des_robot(stateCommand->data.stateDes[6], stateCommand->data.stateDes[7],0);
    Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
    //float trajInitial[12] = {0,0,0, 0,0,.25, 0,0,0,0,0,0};

    if(current_gait == 4)
    {
      float trajInitial[12] = {(float)stateCommand->data.stateDes[3],
                               (float)stateCommand->data.stateDes[4] /*-hw_i->state_estimator->se_ground_pitch*/,
                               (float)stand_traj[5]/*+(float)stateCommand->data.stateDes[11]*/,
                               (float)stand_traj[0]/*+(float)fsm->main_control_settings.p_des[0]*/,
                               (float)stand_traj[1]/*+(float)fsm->main_control_settings.p_des[1]*/,
                               (float)0.26/*fsm->main_control_settings.p_des[2]*/,
                               0,0,0,0,0,0};

      for(int i = 0; i < horizonLength; i++)
        for(int j = 0; j < 12; j++)
          trajAll[12*i+j] = trajInitial[j];
    }

    else
    {
      const float max_pos_error = .1;
      float xStart = world_position_desired[0];
      float yStart = world_position_desired[1];

      //printf("orig \t%.3f\t%.3f\n", xStart, yStart);
      //printf("ref: \t%.3f\t%.3f\n", p[0], p[1]);

      if(xStart - p[0] > max_pos_error) xStart = p[0] + max_pos_error;
      if(p[0] - xStart > max_pos_error) xStart = p[0] - max_pos_error;

      if(yStart - p[1] > max_pos_error) yStart = p[1] + max_pos_error;
      if(p[1] - yStart > max_pos_error) yStart = p[1] - max_pos_error;

      world_position_desired[0] = xStart;
      world_position_desired[1] = yStart;

      //printf("xys: \t%.3f\t%3.f\n", xStart, yStart);
      //std::cout << "pdes_world: " << world_position_desired.transpose() << "\n";
      //printf("perr \t%.3f\t%.3f\n", p[0] - world_position_desired[0], p[1] - world_position_desired[1]);

      float trajInitial[12] = {(float)rpy_comp[0],  // 0
                               (float)rpy_comp[1]/*-hw_i->state_estimator->se_ground_pitch*/,    // 1
                               (float)stateCommand->data.stateDes[5],    // 2
                               xStart,                                   // 3
                               yStart,                                   // 4
                               (float)0.26,      // 5
                               0,                                        // 6
                               0,                                        // 7
                               (float)stateCommand->data.stateDes[11],  // 8
                               v_des_world[0],                           // 9
                               v_des_world[1],                           // 10
                               0};                                       // 11
      for(int i = 0; i < horizonLength; i++)
      {
        for(int j = 0; j < 12; j++)
          trajAll[12*i+j] = trajInitial[j];

        if(i == 0) // start at current position  TODO consider not doing this
        {
          //trajAll[3] = hw_i->state_estimator->se_pBody[0];
          //trajAll[4] = hw_i->state_estimator->se_pBody[1];
          trajAll[2] = seResult.rpy[2];
        }
        else
        {
          trajAll[12*i + 3] = trajAll[12 * (i - 1) + 3] + dtMPC * v_des_world[0];
          trajAll[12*i + 4] = trajAll[12 * (i - 1) + 4] + dtMPC * v_des_world[1];
          trajAll[12*i + 2] = trajAll[12 * (i - 1) + 2] + dtMPC * stateCommand->data.stateDes[11];
        }
      }
    }

//        for(int i = 0; i < 12; i++)
//            printf("%.4f, ", trajAll[i]);

//        printf("\n\n");



    //Timer t1;
    dtMPC = .001 * iterationsBetweenMPC;
    setup_problem(dtMPC,horizonLength,0.4,120);
    //t1.stopPrint("Setup MPC");

    Timer t2;
    //cout << "dtMPC: " << dtMPC << "\n";
    update_problem_data_floats(p,v,q,w,r,yaw,weights,trajAll,alpha,mpcTable);
    //t2.stopPrint("Run MPC");
    //printf("MPC Solve time %f ms\n", t2.getMs());

    for(int leg = 0; leg < 4; leg++)
    {
      Vec3<float> f;
      for(int axis = 0; axis < 3; axis++)
        f[axis] = get_solution(leg*3 + axis);

      f_ff[leg] = -seResult.rBody * f;
    }
  }

}
