#include "FrontJumpCtrl.hpp"


template <typename T>
FrontJumpCtrl<T>::FrontJumpCtrl(DataReader* data_reader,float _dt) : DataReadCtrl<T>(data_reader, _dt) {}


template <typename T>
FrontJumpCtrl<T>::~FrontJumpCtrl() {}

template <typename T>
void FrontJumpCtrl<T>::OneStep(float _curr_time, bool b_preparation, LegControllerCommand<T>* command) {
  DataCtrl::_state_machine_time = _curr_time - DataCtrl::_ctrl_start_time;

  DataCtrl::_b_Preparation = b_preparation;
  _update_joint_command();

  for (int leg = 0; leg < 4; ++leg) {
    for (int jidx = 0; jidx < 3; ++jidx) {
      command[leg].tauFeedForward[jidx] = DataCtrl::_jtorque[3 * leg + jidx];
      command[leg].qDes[jidx] = DataCtrl::_des_jpos[3 * leg + jidx] + 0 * _curr_time;
      command[leg].qdDes[jidx] = DataCtrl::_des_jvel[3 * leg + jidx];
      command[leg].kpJoint(jidx, jidx) = DataCtrl::_Kp_joint[jidx];
      command[leg].kdJoint(jidx, jidx) = DataCtrl::_Kd_joint[jidx];
    }
  }
}

template <typename T>
void FrontJumpCtrl<T>::_update_joint_command() {
  int pre_mode_duration(700);
  int leg_clearance_iteration_front(240) ; 
  //int leg_clearance_iteration_front(180) ; 
  int leg_clearance_iteration(600);
  int leg_ramp_iteration(610);
  int tuck_iteration(610);
  int ramp_end_iteration(700);

  float tau_mult;

  DataCtrl::_des_jpos.setZero();
  DataCtrl::_des_jvel.setZero();
  DataCtrl::_jtorque.setZero();

  DataCtrl::_Kp_joint={10.0, 10.0, 10.0};
  DataCtrl::_Kd_joint={1.0, 1.0, 1.0};

  //DataCtrl::_Kp_joint={20.0, 20.0, 20.0};
  //DataCtrl::_Kd_joint={3.0, 3.0, 3.0};


  // PRE JUMP PREPATATION - CROUCH (FOLLOWS PREMODE DURATION TIMESTEPS) 
  if ( (DataCtrl::pre_mode_count <  pre_mode_duration) || DataCtrl::_b_Preparation) {  
    // move to the initial configuration to prepare for
    // FrontJumping
    if (DataCtrl::pre_mode_count == 0) {
      printf("plan_timesteps: %d \n", DataCtrl::_data_reader->plan_timesteps);
    }
    // printf("pre_mode_count: %d \n", pre_mode_count);
    
    DataCtrl::pre_mode_count += DataCtrl::_key_pt_step;
    DataCtrl::current_iteration = 0;
    tau_mult = 0;
    //DataCtrl::_Kp_joint={20.0, 20.0, 20.0};
    //DataCtrl::_Kd_joint={2.0, 2.0, 2.0};

  } else {
    tau_mult = 1.2;
    // tau_mult = 1.;
  }

  
  // OBTAIN TIMSTEP DATA FROM THE DATA FILE 
  if (DataCtrl::current_iteration > DataCtrl::_data_reader->plan_timesteps - 1) {
    DataCtrl::current_iteration = DataCtrl::_data_reader->plan_timesteps - 1;
  }

  // OBTAIN DATA FROM THE JUMP_DATA FILE GENERATED IN MATLAB 
  float* current_step = DataCtrl::_data_reader->get_plan_at_time(DataCtrl::current_iteration);
  float* tau = current_step + tau_offset;

  // INITIALIZE JOINT PARAMETERS AND TORQUES 
  Vec3<float> q_des_front;
  Vec3<float> q_des_rear;
  Vec3<float> qd_des_front;
  Vec3<float> qd_des_rear;
  Vec3<float> tau_front;
  Vec3<float> tau_rear;

  // SETTING THE JOINT POSITIONS AND VELOCITIES AND FEEDING FORWARD THE JOINT TORQUES obtained from the data file 
  q_des_front << 0.0, current_step[3], current_step[4];
  q_des_rear << 0.0, current_step[5], current_step[6];
  qd_des_front << 0.0, current_step[10], current_step[11];
  qd_des_rear << 0.0, current_step[12], current_step[13];
  tau_front << 0.0, tau_mult * tau[0] / 2.0, tau_mult * tau[1] / 2.0;
  tau_rear << 0.0, tau_mult * tau[2] / 2.0, tau_mult * tau[3] / 2.0;

  // Limitation set so that the arms do not swing too far back and hit the legs 
  if(q_des_front[1] < -M_PI/2.2){
    q_des_front[1] = -M_PI/2.2;
    qd_des_front[1] = 0.;
    tau_front[1] = 0.;
  }
  //pretty_print(tau_front, std::cout, "tau front");
  //pretty_print(tau_rear, std::cout, "tau rear");
  float s(0.);

  
  // CONTROL LEG_CLEARANCE_ITERATION_FRONT
  if (DataCtrl::current_iteration >= leg_clearance_iteration_front &&
      DataCtrl::current_iteration <=leg_clearance_iteration){
    q_des_front << 0.0, current_step[3], current_step[4];

    current_step = DataCtrl::_data_reader->get_plan_at_time(leg_clearance_iteration_front);
    q_des_front << 0.0, -2.3, 2.5;
    qd_des_front << 0.0, 0.0, 0.0;
    tau_front << 0.0, 0.0, 0.0;
  // implement the desried joint states to keep the hands in 
  }
  
  // implement another controller to keep the legs in
  
  // CONTROL LEG_CLEARNACE ITERATION 
  if (DataCtrl::current_iteration >= leg_clearance_iteration 
      && DataCtrl::current_iteration < tuck_iteration) {  // ramp to leg clearance for obstacle
    qd_des_front << 0.0, 0.0, 0.0;
    qd_des_rear << 0.0, 0.0, 0.0;
    tau_front << 0.0, 0.0, 0.0;
    tau_rear << 0.0, 0.0, 0.0;

    s = (float)(DataCtrl::current_iteration - leg_clearance_iteration) /
        (leg_ramp_iteration - leg_clearance_iteration);

    if (s > 1) {
      s = 1;
    }

    Vec3<float> q_des_front_0;
    Vec3<float> q_des_rear_0;
    Vec3<float> q_des_front_f;
    Vec3<float> q_des_rear_f;

    current_step = DataCtrl::_data_reader->get_plan_at_time(tuck_iteration);
    //q_des_front_0 << 0.0, current_step[3], current_step[4];
    q_des_front_0 << 0.0, -2.3, 2.5;
    q_des_rear_0 << 0.0, current_step[5], current_step[6];
    
    // SET THE DESIRED JOINT STATES FOR LEG_CLEARANCE_ITERATION
    current_step = DataCtrl::_data_reader->get_plan_at_time(0);
    q_des_front_f << 0.0, -2.3, 2.5;
    q_des_rear_f << 0.0, -1.25, 2.5;

    // linear interpolation for the ramp 
    q_des_front = (1 - s) * q_des_front_0 + s * q_des_front_f;
    q_des_rear = (1 - s) * q_des_rear_0 + s * q_des_rear_f;
    

  } else if (DataCtrl::current_iteration >= tuck_iteration) { // ramp to landing configuration
    qd_des_front << 0.0, 0.0, 0.0;
    qd_des_rear << 0.0, 0.0, 0.0;
    tau_front << 0.0, 0.0, 0.0;
    tau_rear << 0.0, 0.0, 0.0;

    s = (float)(DataCtrl::current_iteration - tuck_iteration) /
        (ramp_end_iteration - tuck_iteration);

    if (s > 1) {
      s = 1;
    }
    Vec3<float> q_des_front_0;
    Vec3<float> q_des_rear_0;
    Vec3<float> q_des_front_f;
    Vec3<float> q_des_rear_f;

    current_step = DataCtrl::_data_reader->get_plan_at_time(tuck_iteration);
    //q_des_front_0 << 0.0, current_step[3], current_step[4];
    //q_des_rear_0 << 0.0, current_step[5], current_step[6];
    q_des_front_0 << 0.0, -2.3, 2.5;
    q_des_rear_0 << 0.0, -1.25, 2.5;

    // SET THE DESIRED JOINT STATES FOR THE TUCK ITERATION CONTROLLER (LANDING CONTROLLER)
    current_step = DataCtrl::_data_reader->get_plan_at_time(0);
    // q_des_front_f << 0.0, current_step[3], current_step[4];
    // q_des_rear_f << 0.0, current_step[5], current_step[6];
    //q_des_front_f << 0.0, -0.9, 1.8;
    //q_des_front_f << 0.0, -1.0, 2.05;
    //q_des_front_f << 0.0, -0.85, 1.9;
    q_des_front_f << 0.0, -0.45, 1.3;
    
    //q_des_rear_f << 0.0, -0.8, 1.2;
    //q_des_rear_f << 0.0, -0.8, 1.6;
    //q_des_rear_f << 0.0, -1.0, 2.05;
    q_des_rear_f << 0.0, -0.45, 1.3;
    //q_des_rear_f << 0.0, -0.9, 1.8;
   //q_des_front_f << 0.0, -1.2, 2.4;
    //q_des_rear_f << 0.0, -1.2, 2.4;

    // ramp for linear interpolation for the tuck iteration 
    q_des_front = (1 - s) * q_des_front_0 + s * q_des_front_f;
    q_des_rear = (1 - s) * q_des_rear_0 + s * q_des_rear_f;

    DataCtrl::_Kp_joint={20.0, 20.0, 20.0};
    DataCtrl::_Kd_joint={2.0, 2.0, 2.0};
  }

  
  // Abduction
  for (int i = 0; i < 12; i += 3) {
    DataCtrl::_des_jpos[i] = 0.0;
    DataCtrl::_des_jvel[i] = 0.0;
    DataCtrl::_jtorque[i] = 0.0;
  }
  DataCtrl::_des_jpos[0] = s * (-0.2);
  DataCtrl::_des_jpos[3] = s * (0.2);
  DataCtrl::_des_jpos[6] = s * (-0.2);
  DataCtrl::_des_jpos[9] = s * (0.2);

  if(DataCtrl::current_iteration >= tuck_iteration){
  DataCtrl::_des_jpos[0] = (-0.2);
  DataCtrl::_des_jpos[3] = (0.2);
  DataCtrl::_des_jpos[6] = (-0.2);
  DataCtrl::_des_jpos[9] = (0.2);

   }

  // Front Hip
  for (int i = 1; i < 6; i += 3) {
    DataCtrl::_des_jpos[i] = q_des_front[1];
    DataCtrl::_des_jvel[i] = qd_des_front[1];
    DataCtrl::_jtorque[i] = tau_front[1];
  }

  // Front Knee
  for (int i = 2; i < 6; i += 3) {
    DataCtrl::_des_jpos[i] = q_des_front[2];
    DataCtrl::_des_jvel[i] = qd_des_front[2];
    DataCtrl::_jtorque[i] = tau_front[2];
  }

  // Hind Hip
  for (int i = 7; i < 12; i += 3) {
    DataCtrl::_des_jpos[i] = q_des_rear[1];
    DataCtrl::_des_jvel[i] = qd_des_rear[1];
    DataCtrl::_jtorque[i] = tau_rear[1];
  }

  // Hind Knee
  for (int i = 8; i < 12; i += 3) {
    DataCtrl::_des_jpos[i] = q_des_rear[2];
    DataCtrl::_des_jvel[i] = qd_des_rear[2];
    DataCtrl::_jtorque[i] = tau_rear[2];
  }

  // Update rate 0.5kHz
  DataCtrl::current_iteration += DataCtrl::_key_pt_step;
}


template class FrontJumpCtrl<double>;
template class FrontJumpCtrl<float>;
