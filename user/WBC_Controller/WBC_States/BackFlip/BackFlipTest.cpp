#include "BackFlipTest.hpp"

#include <WBC_States/BackFlip/CtrlSet/BackFlipCtrl.hpp>
#include <WBC_States/common/CtrlSet/BodyCtrl.hpp>
#include <WBC_States/common/CtrlSet/FullContactTransCtrl.hpp>
#include <WBC_States/common/CtrlSet/JPosCtrl.hpp>

#include <Math/orientation_tools.h>
#include <Utilities/save_file.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <ParamHandler/ParamHandler.hpp>

template <typename T>
BackFlipTest<T>::BackFlipTest(FloatingBaseModel<T>* robot,
                              const RobotType& type)
    : Test<T>(robot, type) {
  Test<T>::_phase = BFlipPhase::BFlip_body_up_ctrl;
  Test<T>::_state_list.clear();

  _data_reader = new DataReader(type);

  body_up_ctrl_ = new FullContactTransCtrl<T>(robot);
  // backflip_pre_ = new BackFlipCtrl<T>(robot, _data_reader);
  // //JPosCtrl<T>(robot);
  backflip_ctrl_ = new BackFlipCtrl<T>(robot, _data_reader);
  backflip_landing_ =
      new FullContactTransCtrl<T>(robot);  // new JPosCtrl<T>(robot);
  body_ctrl_ = new BodyCtrl<T>(robot);

  Test<T>::_state_list.push_back(body_up_ctrl_);
  // Test<T>::_state_list.push_back(backflip_pre_);
  Test<T>::_state_list.push_back(backflip_ctrl_);
  Test<T>::_state_list.push_back(backflip_landing_);
  Test<T>::_state_list.push_back(body_ctrl_);

  _sp = StateProvider<T>::getStateProvider();

  _SettingParameter();
  if (Test<T>::_b_save_file) {
    _folder_name = "/user/WBC_Controller/sim_data/";
    create_folder(_folder_name);
  }
  printf("[Backflip Test] Constructed\n");
}

template <typename T>
BackFlipTest<T>::~BackFlipTest() {
  for (size_t i(0); i < Test<T>::_state_list.size(); ++i) {
    delete Test<T>::_state_list[i];
  }
}

template <typename T>
void BackFlipTest<T>::_TestInitialization() {
  // Yaml file name
  body_up_ctrl_->CtrlInitialization("CTRL_move_to_target_height");
  // backflip_pre_->CtrlInitialization("CTRL_backflip_pre");
  backflip_ctrl_->CtrlInitialization("CTRL_backflip");
  backflip_landing_->CtrlInitialization("CTRL_backflip_landing");
  body_ctrl_->CtrlInitialization("CTRL_fix_stance");
}

template <typename T>
int BackFlipTest<T>::_NextPhase(const int& phase) {
  int next_phase = phase + 1;
  printf("next phase: %i\n", next_phase);
  if (next_phase == NUM_BFlip_PHASE) {
    return BFlipPhase::BFlip_backflip_end;
  } else
    return next_phase;
}

template <typename T>
void BackFlipTest<T>::_SettingParameter() {
  typename std::vector<Controller<T>*>::iterator iter =
      Test<T>::_state_list.begin();

  while (iter < Test<T>::_state_list.end()) {
    if (Test<T>::_robot_type == RobotType::CHEETAH_3) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_backflip_cheetah3.yaml");
      ParamHandler handler(CheetahConfigPath "TEST_backflip_cheetah3.yaml");
      handler.getBoolean("save_file", Test<T>::_b_save_file);

    } else if (Test<T>::_robot_type == RobotType::MINI_CHEETAH) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_backflip_mini_cheetah.yaml");
      ParamHandler handler(CheetahConfigPath "TEST_backflip_mini_cheetah.yaml");
      handler.getBoolean("save_file", Test<T>::_b_save_file);

    } else {
      printf("[Body Ctrl Test] Invalid robot type\n");
    }
    ++iter;
  }
}

template <typename T>
void BackFlipTest<T>::_UpdateTestOneStep() {
  // Data Save
  if (Test<T>::_b_save_file) {
    static int count(0);
    saveValue(_sp->_curr_time, _folder_name, "time");
    saveVector(_sp->_Q, _folder_name, "config");
    saveVector(_sp->_Qdot, _folder_name, "qdot");
    saveValue(Test<T>::_phase, _folder_name, "phase");
    ++count;
  }
}

template <typename T>
void BackFlipTest<T>::_UpdateExtraData(Cheetah_Extra_Data<T>* ext_data) {
  ext_data->num_step = 0;
  ext_data->num_path_pt = 0;
}

template class BackFlipTest<double>;
template class BackFlipTest<float>;
