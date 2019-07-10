#include "WBICTrotTest.hpp"

#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/CtrlSet/FullContactTransCtrl.hpp>

#include <WBC_States/WBICTrot/CtrlSet/WBIC_FullContactCtrl.hpp>
#include <WBC_States/WBICTrot/CtrlSet/WBIC_TwoContactTransCtrl.hpp>
#include <WBC_States/WBICTrot/CtrlSet/WBIC_TwoLegSwingCtrl.hpp>

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <ParamHandler/ParamHandler.hpp>
#include "Path.hpp"

#include <Utilities/Utilities_print.h>
#include <Utilities/save_file.h>

template <typename T>
WBICTrotTest<T>::WBICTrotTest(FloatingBaseModel<T> *robot,
                              const RobotType &type)
    : Test<T>(robot, type),
      _des_jpos(cheetah::num_act_joint),
      _des_jvel(cheetah::num_act_joint),
      _des_jacc(cheetah::num_act_joint),
      _b_remote_ctrl(true),
      _visionLCM(getLcmUrl(255)),
      _b_loc_update(false),
      _wbcLCM(getLcmUrl(255)) {
  _vision_loc.setZero();
  _visionLCM.subscribe("vision_data", &WBICTrotTest<T>::handleVisionLCM, this);
  _visionLCMThread = std::thread(&WBICTrotTest<T>::visionLCMThread, this);

  //_path = new LinPath<T>();
  _path = new CircularPath<T>();

  _body_pos.setZero();
  _body_vel.setZero();
  _body_acc.setZero();

  _body_ori_rpy.setZero();
  _body_ang_vel.setZero();

  Test<T>::_phase = WBICTrotPhase::lift_up;
  Test<T>::_state_list.clear();

  body_up_ctrl_ = new FullContactTransCtrl<T>(robot);
  body_ctrl_ = new WBIC_FullContactCtrl<T>(this, robot);

  frhl_swing_start_trans_ctrl_ =
      new WBIC_TwoContactTransCtrl<T>(this, robot, linkID::FR, linkID::HL, 1);
  frhl_swing_ctrl_ =
      new WBIC_TwoLegSwingCtrl<T>(this, robot, linkID::FR, linkID::HL);
  frhl_swing_end_trans_ctrl_ =
      new WBIC_TwoContactTransCtrl<T>(this, robot, linkID::FR, linkID::HL, -1);

  flhr_swing_start_trans_ctrl_ =
      new WBIC_TwoContactTransCtrl<T>(this, robot, linkID::FL, linkID::HR, 1);
  flhr_swing_ctrl_ =
      new WBIC_TwoLegSwingCtrl<T>(this, robot, linkID::FL, linkID::HR);
  flhr_swing_end_trans_ctrl_ =
      new WBIC_TwoContactTransCtrl<T>(this, robot, linkID::FL, linkID::HR, -1);

  Test<T>::_state_list.push_back(body_up_ctrl_);
  Test<T>::_state_list.push_back(body_ctrl_);

  Test<T>::_state_list.push_back(frhl_swing_start_trans_ctrl_);
  Test<T>::_state_list.push_back(frhl_swing_ctrl_);
  Test<T>::_state_list.push_back(frhl_swing_end_trans_ctrl_);

  Test<T>::_state_list.push_back(body_ctrl_);

  Test<T>::_state_list.push_back(flhr_swing_start_trans_ctrl_);
  Test<T>::_state_list.push_back(flhr_swing_ctrl_);
  Test<T>::_state_list.push_back(flhr_swing_end_trans_ctrl_);

  _sp = StateProvider<T>::getStateProvider();
  _SettingParameter();

  // (x, y)
  _filtered_input_vel.push_back(
      new digital_lp_filter<T>(2. * M_PI * 15., Test<T>::dt));
  _filtered_input_vel.push_back(
      new digital_lp_filter<T>(2. * M_PI * 15., Test<T>::dt));

  _input_vel.setZero();

  if (Test<T>::_b_save_file) {
    _folder_name = "/user/WBC_Controller/sim_data/";
    create_folder(_folder_name);
  }
  printf("[WBIC Trot Test] Constructed\n");
}

template <typename T>
void WBICTrotTest<T>::handleVisionLCM(const lcm::ReceiveBuffer *rbuf,
                                      const std::string &chan,
                                      const vision_data_t *msg) {
  (void)rbuf;
  (void)chan;
  printf("lidar pos: %f, %f, %f\n", msg->p_lidar[0], msg->p_lidar[1],
         msg->p_lidar[2]);
  printf("lidar quaternion (x, y, z, w) : %f, %f, %f, %f\n",
         msg->lidar_quaternion[0], msg->lidar_quaternion[1],
         msg->lidar_quaternion[2], msg->lidar_quaternion[3]);

  for (size_t i(0); i < 2; ++i) _vision_loc[i] = msg->p_lidar[i];
  Quat<T> vision_quat;
  vision_quat[0] = msg->lidar_quaternion[3];
  vision_quat[1] = msg->lidar_quaternion[0];
  vision_quat[2] = msg->lidar_quaternion[1];
  vision_quat[3] = msg->lidar_quaternion[2];

  Vec3<T> rpy = ori::quatToRPY(vision_quat);
  _vision_loc[2] = rpy[2];
  //_vision_loc[0] = -msg->p_lidar[2];
  //_vision_loc[1] = -msg->p_lidar[0];
  //_vision_loc[2] = msg->p_lidar[1];
}

template <typename T>
void WBICTrotTest<T>::visionLCMThread() {
  while (true) {
    _visionLCM.handle();
    _b_loc_update = true;
  }
}

template <typename T>
WBICTrotTest<T>::~WBICTrotTest() {
  for (size_t i(0); i < Test<T>::_state_list.size(); ++i) {
    delete Test<T>::_state_list[i];
  }
}

template <typename T>
void WBICTrotTest<T>::_TestInitialization() {
  // Yaml file name
  body_up_ctrl_->CtrlInitialization("CTRL_move_to_target_height");
  body_ctrl_->CtrlInitialization("CTRL_fix_stance");

  // Transition
  frhl_swing_start_trans_ctrl_->CtrlInitialization("CTRL_two_leg_trans");
  frhl_swing_end_trans_ctrl_->CtrlInitialization("CTRL_two_leg_trans");
  flhr_swing_start_trans_ctrl_->CtrlInitialization("CTRL_two_leg_trans");
  flhr_swing_end_trans_ctrl_->CtrlInitialization("CTRL_two_leg_trans");

  // Swing
  frhl_swing_ctrl_->CtrlInitialization("CTRL_frhl_swing");
  flhr_swing_ctrl_->CtrlInitialization("CTRL_flhr_swing");
}

template <typename T>
int WBICTrotTest<T>::_NextPhase(const int &phase) {
  int next_phase = phase + 1;
  if (next_phase == WBICTrotPhase::NUM_TROT_PHASE) {
    next_phase = WBICTrotPhase::full_contact_1;
  }
  // printf("next phase: %i\n", next_phase);

  // Which path segment we are going to use in this swing
  if (next_phase == WBICTrotPhase::flhr_swing_start_trans ||
      next_phase == WBICTrotPhase::frhl_swing_start_trans) {
  }
  // Check whether Cheetah must stop or step
  if (next_phase == WBICTrotPhase::flhr_swing_start_trans) {
    // pretty_print(_front_foot_loc, std::cout, "front foot");
    // pretty_print(_hind_foot_loc, std::cout, "hind foot");
    Vec3<T> landing_loc_ave = Vec3<T>::Zero();
    landing_loc_ave += 0.5 * Test<T>::_robot->_pGC[linkID::FR];
    landing_loc_ave += 0.5 * Test<T>::_robot->_pGC[linkID::HL];

    landing_loc_ave[2] = 0.;  // Flat terrain

    if (!_b_loc_update) {
      _vision_loc += landing_loc_ave;
      Vec3<T> rpy = ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      _vision_loc[2] = rpy[2];
    }

    _body_pos.setZero();
    _body_pos[2] = _target_body_height;

    _sp->_contact_pt[0] = linkID::FR;
    _sp->_contact_pt[1] = linkID::HL;
    _sp->_num_contact = 2;

    //_sp->_local_frame_global_pos = landing_loc_ave;
    _sp->_local_frame_global_pos.setZero();
    // pretty_print(_sp->_local_frame_global_pos, std::cout, "local frame");
  }

  if (next_phase == WBICTrotPhase::frhl_swing_start_trans) {
    Vec3<T> landing_loc_ave = Vec3<T>::Zero();
    landing_loc_ave += 0.5 * Test<T>::_robot->_pGC[linkID::FL];
    landing_loc_ave += 0.5 * Test<T>::_robot->_pGC[linkID::HR];

    landing_loc_ave[2] = 0.;  // Flat terrain

    if (!_b_loc_update) {
      _vision_loc += landing_loc_ave;
      Vec3<T> rpy = ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      _vision_loc[2] = rpy[2];
    }

    _body_pos.setZero();
    _body_pos[2] = _target_body_height;

    _sp->_contact_pt[0] = linkID::FL;
    _sp->_contact_pt[1] = linkID::HR;
    _sp->_num_contact = 2;

    //_sp->_local_frame_global_pos = landing_loc_ave;
    _sp->_local_frame_global_pos.setZero();
    // pretty_print(_sp->_local_frame_global_pos, std::cout, "local frame");
  }
  return next_phase;
}

template <typename T>
void WBICTrotTest<T>::_SettingParameter() {
  typename std::vector<Controller<T> *>::iterator iter =
      Test<T>::_state_list.begin();
  ParamHandler *handler = NULL;

  if (Test<T>::_robot_type == RobotType::CHEETAH_3) {
    while (iter < Test<T>::_state_list.end()) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_wbic_trot_cheetah3.yaml");
      ++iter;
    }
    handler =
        new ParamHandler(CheetahConfigPath "TEST_wbic_trot_cheetah3.yaml");
  } else if (Test<T>::_robot_type == RobotType::MINI_CHEETAH) {
    while (iter < Test<T>::_state_list.end()) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_wbic_trot_mini_cheetah.yaml");
      ++iter;
    }
    handler =
        new ParamHandler(CheetahConfigPath "TEST_wbic_trot_mini_cheetah.yaml");
  } else {
    printf("[WBIC Trot Test] Invalid robot type\n");
  }

  handler->getValue<T>("body_height", _target_body_height);
  _body_pos[2] = _target_body_height;
  handler->getBoolean("remote_control", _b_remote_ctrl);
  handler->getBoolean("save_file", Test<T>::_b_save_file);
  delete handler;
}

template <typename T>
void WBICTrotTest<T>::_UpdateTestOneStep() {
  static int count(0);
  ++count;
  if (_b_remote_ctrl) {
    _body_ang_vel[2] = _sp->_ori_command[2];

    // Command Truncation
    if (fabs(_body_ang_vel[2]) < 0.1) {
      _body_ang_vel[2] = 0.;
    }
    if (fabs(_sp->_dir_command[0]) < 0.1) {
      _sp->_dir_command[0] = 0.;
    }
    if (fabs(_sp->_dir_command[1]) < 0.1) {
      _sp->_dir_command[1] = 0.;
    }

    _filtered_input_vel[0]->input(_sp->_dir_command[0]);
    _filtered_input_vel[1]->input(_sp->_dir_command[1]);

    _input_vel[0] = _filtered_input_vel[0]->output();
    _input_vel[1] = _filtered_input_vel[1]->output();

    _body_ori_rpy[2] += _body_ang_vel[2] * Test<T>::dt;
    Mat3<T> Rot = rpyToRotMat(_body_ori_rpy);
    _body_vel = Rot.transpose() * _input_vel;
    _body_pos += _body_vel * Test<T>::dt;

  } else {  // Vision Data
    Vec3<T> des_loc;
    des_loc.setZero();
    if (_sp->_mode == 11) {
      _path->getDesLoc(_sp->_curr_time, _vision_loc, des_loc);

      _input_vel[0] = 1.3 * (des_loc[0] - _vision_loc[0]);
      _input_vel[1] = 1.3 * (des_loc[1] - _vision_loc[1]);
      // T yaw_err =(des_loc[2] - _vision_loc[2]);
      // TEST
      Vec3<T> rpy = ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      T yaw_err = (des_loc[2] - rpy[2]);

      while (true) {
        if (yaw_err > M_PI) {
          yaw_err -= 2. * M_PI;
        } else if (yaw_err < -M_PI) {
          yaw_err += 2. * M_PI;
        } else {
          break;
        }
      }
      _body_ang_vel[2] = 2. * yaw_err;
    }
    // printf("yaw cmd, act, err: %f, %f, %f\n", des_loc[2],
    //_vision_loc[2], yaw_err);

    _body_ang_vel[2] = coerce<T>(_body_ang_vel[2], -0.6, 0.6);
    _input_vel[0] = coerce<T>(_input_vel[0], -0.5, 0.5);
    _input_vel[1] = coerce<T>(_input_vel[1], -0.5, 0.5);

    _body_ori_rpy[2] += _body_ang_vel[2] * Test<T>::dt;
    _body_vel = _input_vel;
    //_body_pos += _body_vel*Test<T>::dt;
    if (count % 300 == 0) {
      printf("target loc: %f, %f\n", des_loc[0], des_loc[1]);
      printf("vision loc: %f, %f\n", _vision_loc[0], _vision_loc[1]);
      printf("loc err: %f, %f\n", des_loc[0] - _vision_loc[0],
             des_loc[1] - _vision_loc[1]);
    }
  }

  // printf("command: %f, %f, %f \n", _input_vel[0], _input_vel[1],
  // _body_ang_vel[2]);

  // T body_height_adj = 0.05 * sqrt(_body_vel[0]*_body_vel[0] +
  // _body_vel[1]*_body_vel[1]); body_height_adj = coerce<T>(body_height_adj, 0.,
  // 0.2); _body_pos[2] = _target_body_height + body_height_adj; _body_pos =
  //Test<T>::_robot->_state.bodyPosition;
}

template <typename T>
void WBICTrotTest<T>::_UpdateExtraData(Cheetah_Extra_Data<T> *ext_data) {
  (void)ext_data;

  if (Test<T>::_b_save_file) {
    static int count(0);
    if (count % 10 == 0) {
      saveValue(_sp->_curr_time, _folder_name, "time");
      saveVector(_body_pos, _folder_name, "body_pos");
      saveVector(_body_vel, _folder_name, "body_vel");
      saveVector(_body_acc, _folder_name, "body_acc");

      Vec3<T> body_ori_rpy =
          ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      saveVector(body_ori_rpy, _folder_name, "body_ori_rpy");
      saveVector(_body_ang_vel, _folder_name, "body_ang_vel");
      saveVector(_body_ori_rpy, _folder_name, "cmd_body_ori_rpy");

      saveVector(_sp->_Q, _folder_name, "config");
      saveVector(_sp->_Qdot, _folder_name, "qdot");
      saveVector((Test<T>::_copy_cmd)[0].tauFeedForward, _folder_name,
                 "fr_tau");
      saveVector((Test<T>::_copy_cmd)[1].tauFeedForward, _folder_name,
                 "fl_tau");
      saveVector((Test<T>::_copy_cmd)[2].tauFeedForward, _folder_name,
                 "hr_tau");
      saveVector((Test<T>::_copy_cmd)[3].tauFeedForward, _folder_name,
                 "hl_tau");

      saveValue(Test<T>::_phase, _folder_name, "phase");
    }
    ++count;
  }

  // LCM
  for (size_t i(0); i < cheetah::num_act_joint; ++i) {
    _wbc_data_lcm.jpos_cmd[i] = _des_jpos[i];
    _wbc_data_lcm.jvel_cmd[i] = _des_jvel[i];
    _wbc_data_lcm.jacc_cmd[i] = _des_jacc[i];
  }

  //_wbc_data_lcm.body_ori[3] = Ctrl::_robot_sys->_state.bodyOrientation[3];
  for (size_t i(0); i < 3; ++i) {
    _wbc_data_lcm.vision_loc[i] = _vision_loc[i];
  }
  //_wbc_data_lcm.body_ori[i] = Ctrl::_robot_sys->_state.bodyOrientation[i];
  //_wbc_data_lcm.body_ang_vel[i] = Ctrl::_robot_sys->_state.bodyVelocity[i];

  //_wbc_data_lcm.fr_foot_pos_cmd[i] = _fr_foot_pos[i];
  //_wbc_data_lcm.fl_foot_pos_cmd[i] = _fl_foot_pos[i];
  //_wbc_data_lcm.hr_foot_pos_cmd[i] = _hr_foot_pos[i];
  //_wbc_data_lcm.hl_foot_pos_cmd[i] = _hl_foot_pos[i];

  //_wbc_data_lcm.fr_foot_vel_cmd[i] = _fr_foot_vel[i];
  //_wbc_data_lcm.fl_foot_vel_cmd[i] = _fl_foot_vel[i];
  //_wbc_data_lcm.hr_foot_vel_cmd[i] = _hr_foot_vel[i];
  //_wbc_data_lcm.hl_foot_vel_cmd[i] = _hl_foot_vel[i];

  //_wbc_data_lcm.fr_foot_local_pos[i] = _fr_foot_pos[i] -
  //_fr_foot_local_task->getPosError()[i]; _wbc_data_lcm.fl_foot_local_pos[i] =
  //_fl_foot_pos[i] - _fl_foot_local_task->getPosError()[i];
  //_wbc_data_lcm.hr_foot_local_pos[i] = _hr_foot_pos[i] -
  //_hr_foot_local_task->getPosError()[i]; _wbc_data_lcm.hl_foot_local_pos[i] =
  //_hl_foot_pos[i] - _hl_foot_local_task->getPosError()[i];

  //_wbc_data_lcm.fr_foot_local_vel[i] =
  //(Ctrl::_robot_sys->_vGC[linkID::FR])[i] -
  //(Ctrl::_robot_sys->_vGC[linkID::FR_abd])[i];
  //_wbc_data_lcm.fl_foot_local_vel[i] =
  //(Ctrl::_robot_sys->_vGC[linkID::FL])[i] -
  //(Ctrl::_robot_sys->_vGC[linkID::FL_abd])[i];
  //_wbc_data_lcm.hr_foot_local_vel[i] =
  //(Ctrl::_robot_sys->_vGC[linkID::HR])[i] -
  //(Ctrl::_robot_sys->_vGC[linkID::HR_abd])[i];
  //_wbc_data_lcm.hl_foot_local_vel[i] =
  //(Ctrl::_robot_sys->_vGC[linkID::HL])[i] -
  //(Ctrl::_robot_sys->_vGC[linkID::HL_abd])[i];

  // if((!_b_front_swing) && (!_b_hind_swing)){
  //_wbc_data_lcm.fr_Fr_des[i] = _wbic_data->_Fr_des[i];
  //_wbc_data_lcm.fl_Fr_des[i] = _wbic_data->_Fr_des[i+3];
  //_wbc_data_lcm.hr_Fr_des[i] = _wbic_data->_Fr_des[i+6];
  //_wbc_data_lcm.hl_Fr_des[i] = _wbic_data->_Fr_des[i+9];

  //_wbc_data_lcm.fr_Fr[i] = _wbic_data->_Fr[i];
  //_wbc_data_lcm.fl_Fr[i] = _wbic_data->_Fr[i+3];
  //_wbc_data_lcm.hr_Fr[i] = _wbic_data->_Fr[i+6];
  //_wbc_data_lcm.hl_Fr[i] = _wbic_data->_Fr[i+9];
  //}
  // else if((!_b_front_swing)){
  //_wbc_data_lcm.fr_Fr_des[i] = _wbic_data->_Fr_des[i];
  //_wbc_data_lcm.fl_Fr_des[i] = _wbic_data->_Fr_des[i+3];

  //_wbc_data_lcm.fr_Fr[i] = _wbic_data->_Fr[i];
  //_wbc_data_lcm.fl_Fr[i] = _wbic_data->_Fr[i+3];
  //}
  // else if((!_b_hind_swing)){
  //_wbc_data_lcm.hr_Fr_des[i] = _wbic_data->_Fr_des[i];
  //_wbc_data_lcm.hl_Fr_des[i] = _wbic_data->_Fr_des[i+3];

  //_wbc_data_lcm.hr_Fr[i] = _wbic_data->_Fr[i];
  //_wbc_data_lcm.hl_Fr[i] = _wbic_data->_Fr[i+3];
  //}
  //}
  _wbcLCM.publish("wbc_lcm_data", &_wbc_data_lcm);
}

template class WBICTrotTest<double>;
template class WBICTrotTest<float>;
