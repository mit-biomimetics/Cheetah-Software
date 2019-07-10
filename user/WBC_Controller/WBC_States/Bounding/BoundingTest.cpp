#include "BoundingTest.hpp"

#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/CtrlSet/FullContactTransCtrl.hpp>
#include <WBC_States/common/CtrlSet/PostureKeepingCtrl.hpp>

#include <WBC_States/Bounding/CtrlSet/KinBoundingCtrl.hpp>
#include <WBC_States/Bounding/CtrlSet/BoundingJumpCtrl.hpp>

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <ParamHandler/ParamHandler.hpp>

#include <Utilities/Utilities_print.h>
#include <Utilities/save_file.h>

template <typename T>
BoundingTest<T>::BoundingTest(FloatingBaseModel<T>* robot,
                              const RobotType& type)
    : Test<T>(robot, type), _stop(false) {
  _body_pos.setZero();
  _body_vel.setZero();
  _body_acc.setZero();

  _body_ori_rpy.setZero();
  _body_ang_vel.setZero();

  Test<T>::_phase = BoundingPhase::lift_up;
  Test<T>::_state_list.clear();

  _body_up_ctrl = new FullContactTransCtrl<T>(robot);
  _posture_keeping = new PostureKeepingCtrl<T>(robot);
  _bounding = new KinBoundingCtrl<T>(this, robot);
  _bounding_jump = new BoundingJumpCtrl<T>(this, robot);

  Test<T>::_state_list.push_back(_body_up_ctrl);
  Test<T>::_state_list.push_back(_posture_keeping);
  Test<T>::_state_list.push_back(_bounding);
  Test<T>::_state_list.push_back(_bounding_jump);

  _sp = StateProvider<T>::getStateProvider();
  _SettingParameter();

  if (Test<T>::_b_save_file) {
    _folder_name = "/user/WBC_Controller/sim_data/";
    create_folder(_folder_name);
  }
  printf("[Bounding Test] Constructed\n");
}

template <typename T>
BoundingTest<T>::~BoundingTest() {
  for (size_t i(0); i < Test<T>::_state_list.size(); ++i) {
    delete Test<T>::_state_list[i];
  }
}

template <typename T>
void BoundingTest<T>::_TestInitialization() {
  // Yaml file name
  _body_up_ctrl->CtrlInitialization("CTRL_move_to_target_height");
  _posture_keeping->CtrlInitialization("CTRL_keep_posture");
  _bounding->CtrlInitialization("CTRL_bounding");
  _bounding_jump->CtrlInitialization("CTRL_bounding_jump");
}

template <typename T>
int BoundingTest<T>::_NextPhase(const int& phase) {
  int next_phase = phase + 1;
  _sp->_num_contact = 0;

  if( _stop ){
    next_phase = BoundingPhase::posture_keeping;
  }else{
    if(phase == BoundingPhase::bounding){
      next_phase = BoundingPhase::bounding_jump;
    }else if(phase == BoundingPhase::bounding_jump){
      next_phase = BoundingPhase::bounding;
    }
  }
  printf("next phase: %i\n", next_phase);
  return next_phase;
}

template <typename T>
void BoundingTest<T>::_SettingParameter() {
  typename std::vector<Controller<T>*>::iterator iter =
      Test<T>::_state_list.begin();
  ParamHandler* handler = NULL;
  while (iter < Test<T>::_state_list.end()) {
    // Cheetah 3
    if (Test<T>::_robot_type == RobotType::CHEETAH_3) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_bounding_cheetah3.yaml");
      handler =
          new ParamHandler(CheetahConfigPath "TEST_bounding_cheetah3.yaml");
      // Mini Cheetah
    } else if (Test<T>::_robot_type == RobotType::MINI_CHEETAH) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_bounding_mini_cheetah.yaml");
      handler =
          new ParamHandler(CheetahConfigPath "TEST_bounding_mini_cheetah.yaml");
    } else {
      printf("[Bounding Test] Invalid robot type\n");
    }
    ++iter;
  }
  handler->getBoolean("save_file", Test<T>::_b_save_file);
  delete handler;
}

template <typename T>
void BoundingTest<T>::_UpdateTestOneStep() {
  // T scale(5.0/3.0);
  T scale(1.0);
  _body_ang_vel[2] = _sp->_ori_command[2];
  _body_ori_rpy[2] += _body_ang_vel[2] * Test<T>::dt;

  Vec3<T> input_vel;
  input_vel.setZero();
  input_vel[0] = scale * _sp->_dir_command[0];
  input_vel[1] = 0.5 * scale * _sp->_dir_command[1];
  // Mat3<T> Rot = rpyToRotMat(_body_ori_rpy);
  //_body_vel = Rot.transpose() * input_vel;
  _body_vel = input_vel;
  //_body_vel[0] = 1.3;
  //_body_pos += _body_vel*Test<T>::dt;
  // printf("input vel: %f, %f\n", input_vel[0], input_vel[1]);
}

template <typename T>
void BoundingTest<T>::_UpdateExtraData(Cheetah_Extra_Data<T>* ext_data) {
  (void)ext_data;

  if (Test<T>::_b_save_file) {
    static int count(0);
    if (count % 10 == 0) {
      saveValue(_sp->_curr_time, _folder_name, "time");
      saveVector(_body_pos, _folder_name, "body_pos");
      saveVector(_body_vel, _folder_name, "body_vel");

      Vec3<T> body_ori_rpy =
          ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      saveVector(body_ori_rpy, _folder_name, "body_ori_rpy");
      saveVector(_body_ori_rpy, _folder_name, "cmd_body_ori_rpy");
      saveVector(_body_ang_vel, _folder_name, "body_ang_vel");

      saveValue(Test<T>::_phase, _folder_name, "phase");

      saveVector(Test<T>::_robot->_pGC[linkID::FR], _folder_name, "fr_pos");
      saveVector(Test<T>::_robot->_pGC[linkID::FL], _folder_name, "fl_pos");
      saveVector(Test<T>::_robot->_pGC[linkID::HR], _folder_name, "hr_pos");
      saveVector(Test<T>::_robot->_pGC[linkID::HL], _folder_name, "hl_pos");

      saveVector(Test<T>::_robot->_vGC[linkID::FR], _folder_name, "fr_vel");
      saveVector(Test<T>::_robot->_vGC[linkID::FL], _folder_name, "fl_vel");
      saveVector(Test<T>::_robot->_vGC[linkID::HR], _folder_name, "hr_vel");
      saveVector(Test<T>::_robot->_vGC[linkID::HL], _folder_name, "hl_vel");

      // saveVector((Test<T>::_copy_cmd)[0].tauFeedForward, _folder_name,
      // "fr_tau"); saveVector((Test<T>::_copy_cmd)[1].tauFeedForward,
      // _folder_name, "fl_tau");
      // saveVector((Test<T>::_copy_cmd)[2].tauFeedForward, _folder_name,
      // "hr_tau"); saveVector((Test<T>::_copy_cmd)[3].tauFeedForward,
      // _folder_name, "hl_tau");
    }
    ++count;
  }
}

template class BoundingTest<double>;
template class BoundingTest<float>;
