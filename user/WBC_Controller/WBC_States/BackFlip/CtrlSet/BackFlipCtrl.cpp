#include "BackFlipCtrl.hpp"

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>
#include <WBC_States/common/TaskSet/BodyOriTask.hpp>
#include <WBC_States/common/TaskSet/BodyPosTask.hpp>
#include <WBC_States/common/TaskSet/LinkPosTask.hpp>

#include <ParamHandler/ParamHandler.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC/WBLC/WBLC.hpp>
#include <WBC_States/Test.hpp>

template <typename T>
BackFlipCtrl<T>::BackFlipCtrl(const FloatingBaseModel<T>* robot,
                              DataReader* data_reader)
    : Controller<T>(robot),
      _data_reader(data_reader),
      _Kp(cheetah::num_act_joint),
      _Kd(cheetah::num_act_joint),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _jtorque(cheetah::num_act_joint),
      _end_time(3.5),
      _dim_contact(0),
      _ctrl_start_time(0.) {
  _sp = StateProvider<T>::getStateProvider();

  printf("[Body Control] Constructed\n");
}

template <typename T>
BackFlipCtrl<T>::~BackFlipCtrl() {
  delete _param_handler;
}

template <typename T>
void BackFlipCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;

  _update_joint_command();

  for (size_t leg(0); leg < cheetah::num_leg; ++leg) {
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
      ((LegControllerCommand<T>*)_cmd)[leg].tauFeedForward[jidx] =
          _jtorque[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qDes[jidx] =
          _des_jpos[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].qdDes[jidx] =
          _des_jvel[cheetah::num_leg_joint * leg + jidx];

      ((LegControllerCommand<T>*)_cmd)[leg].kpJoint(jidx, jidx) =
          _Kp_joint[jidx];
      ((LegControllerCommand<T>*)_cmd)[leg].kdJoint(jidx, jidx) =
          _Kd_joint[jidx];
    }
  }
  Ctrl::_PostProcessing_Command();
}

template <typename T>
void BackFlipCtrl<T>::_update_joint_command() {
  static int current_iteration(0);
  static int pre_mode_count(0);
  int pre_mode_duration(2000);
  int tuck_iteration(600);
  int ramp_end_iteration(650);
  float tau_mult;

  _des_jpos.setZero();
  _des_jvel.setZero();
  _jtorque.setZero();

  if (pre_mode_count <
      pre_mode_duration) {  // move to the initial configuration to prepare for
                            // backfliping
    if (pre_mode_count == 0) {
      printf("plan_timesteps: %d \n", _data_reader->plan_timesteps);
    }
    // printf("pre_mode_count: %d \n", pre_mode_count);
    pre_mode_count++;
    current_iteration = 0;
    tau_mult = 0;
  } else {
    tau_mult = 1.2;
    // tau_mult = 1.;
  }

  if (current_iteration > _data_reader->plan_timesteps - 1) {
    current_iteration = _data_reader->plan_timesteps - 1;
  }

  float* current_step = _data_reader->get_plan_at_time(current_iteration);
  float* tau = current_step + tau_offset;

  Vec3<float> q_des_front;
  Vec3<float> q_des_rear;
  Vec3<float> qd_des_front;
  Vec3<float> qd_des_rear;
  Vec3<float> tau_front;
  Vec3<float> tau_rear;

  q_des_front << 0.0, current_step[3], current_step[4];
  q_des_rear << 0.0, current_step[5], current_step[6];
  qd_des_front << 0.0, current_step[10], current_step[11];
  qd_des_rear << 0.0, current_step[12], current_step[13];
  tau_front << 0.0, tau_mult * tau[0] / 2.0, tau_mult * tau[1] / 2.0;
  tau_rear << 0.0, tau_mult * tau[2] / 2.0, tau_mult * tau[3] / 2.0;

  float s(0.);

  if (current_iteration >= tuck_iteration) {  // ramp to landing configuration
    qd_des_front << 0.0, 0.0, 0.0;
    qd_des_rear << 0.0, 0.0, 0.0;
    tau_front << 0.0, 0.0, 0.0;
    tau_rear << 0.0, 0.0, 0.0;

    s = (float)(current_iteration - tuck_iteration) /
        (ramp_end_iteration - tuck_iteration);

    // printf("tuck ramping phase s: %f \n",s);

    if (s > 1) {
      s = 1;
    }

    Vec3<float> q_des_front_0;
    Vec3<float> q_des_rear_0;
    Vec3<float> q_des_front_f;
    Vec3<float> q_des_rear_f;

    current_step = _data_reader->get_plan_at_time(tuck_iteration);
    q_des_front_0 << 0.0, current_step[3], current_step[4];
    q_des_rear_0 << 0.0, current_step[5], current_step[6];

    current_step = _data_reader->get_plan_at_time(0);
    // q_des_front_f << 0.0, current_step[3], current_step[4];
    // q_des_rear_f << 0.0, current_step[5], current_step[6];
    q_des_front_f << 0.0, -0.8425, 1.65;
    q_des_rear_f << 0.0, -0.8425, 1.65;

    q_des_front = (1 - s) * q_des_front_0 + s * q_des_front_f;
    q_des_rear = (1 - s) * q_des_rear_0 + s * q_des_rear_f;

    // for(int i=0; i<12; i++){
    //	_Kp_joint[i]=10.0;
    //	_Kd_joint[i]=1.0;
    //}
  }

  // Abduction
  for (int i = 0; i < 12; i += 3) {
    _des_jpos[i] = 0.0;
    _des_jvel[i] = 0.0;
    _jtorque[i] = 0.0;
  }
  _des_jpos[0] = s * (-0.2);
  _des_jpos[3] = s * (0.2);
  _des_jpos[6] = s * (-0.2);
  _des_jpos[9] = s * (0.2);

  // Front Hip
  for (int i = 1; i < 6; i += 3) {
    _des_jpos[i] = q_des_front[1];
    _des_jvel[i] = qd_des_front[1];
    _jtorque[i] = tau_front[1];
  }

  // Front Knee
  for (int i = 2; i < 6; i += 3) {
    _des_jpos[i] = q_des_front[2];
    _des_jvel[i] = qd_des_front[2];
    _jtorque[i] = tau_front[2];
  }

  // Hind Hip
  for (int i = 7; i < 12; i += 3) {
    _des_jpos[i] = q_des_rear[1];
    _des_jvel[i] = qd_des_rear[1];
    _jtorque[i] = tau_rear[1];
  }

  // Hind Knee
  for (int i = 8; i < 12; i += 3) {
    _des_jpos[i] = q_des_rear[2];
    _des_jvel[i] = qd_des_rear[2];
    _jtorque[i] = tau_rear[2];
  }

  current_iteration++;
}

template <typename T>
void BackFlipCtrl<T>::FirstVisit() {
  _ctrl_start_time = _sp->_curr_time;
}

template <typename T>
void BackFlipCtrl<T>::LastVisit() {}

template <typename T>
bool BackFlipCtrl<T>::EndOfPhase() {
  if (Ctrl::_state_machine_time > (_end_time - 2. * Test<T>::dt)) {
    return true;
  }

   for(int i=2; i<12; i+=3){ // check contact for each leg

  	 if(Ctrl::_state_machine_time>2.7 && _sp->_Q[i+6]>_q_knee_max && _sp->_Qdot[i+6]>_qdot_knee_max){
  	 	printf("Contact detected at leg [%d] => Switch to the landing phase !!! \n", i/3); 
		printf("state_machine_time: %lf \n",Ctrl::_state_machine_time); 
		printf("Q-Knee: %lf \n",_sp->_Q[i+6]);
  	 	printf("Qdot-Knee: %lf \n",_sp->_Qdot[i+6]);
  	 	return true;
  	}
  }

  return false;
}

template <typename T>
void BackFlipCtrl<T>::CtrlInitialization(const std::string& category_name) {
  _param_handler->getValue<T>(category_name, "qdot_knee_max", _q_knee_max);
  _param_handler->getValue<T>(category_name, "qdot_knee_max", _qdot_knee_max);
}

template <typename T>
void BackFlipCtrl<T>::SetTestParameter(const std::string& test_file) {
  _param_handler = new ParamHandler(test_file);

  std::vector<T> tmp_vec;
  // Feedback Gain
  _param_handler->getVector<T>("Kp", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _Kp[i] = tmp_vec[i];
  }
  _param_handler->getVector<T>("Kd", tmp_vec);
  for (size_t i(0); i < tmp_vec.size(); ++i) {
    _Kd[i] = tmp_vec[i];
  }
  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint_backflip", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint_backflip", _Kd_joint);
}

template class BackFlipCtrl<double>;
template class BackFlipCtrl<float>;
