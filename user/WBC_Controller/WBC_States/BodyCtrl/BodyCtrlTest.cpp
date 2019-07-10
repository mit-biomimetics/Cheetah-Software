#include "BodyCtrlTest.hpp"

#include <Math/orientation_tools.h>
#include <Utilities/save_file.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <ParamHandler/ParamHandler.hpp>
#include <WBC_States/BodyCtrl/BodyPostureCtrl.hpp>
#include <WBC_States/common/CtrlSet/FullContactTransCtrl.hpp>

template <typename T>
BodyCtrlTest<T>::BodyCtrlTest(FloatingBaseModel<T>* robot,
                              const RobotType& type)
    : Test<T>(robot, type) {
  //_phase = BodyCtrlPhase::BDCTRL_body_ctrl;
  Test<T>::_phase = BodyCtrlPhase::BDCTRL_body_up_ctrl;
  Test<T>::_state_list.clear();

  body_up_ctrl_ = new FullContactTransCtrl<T>(robot);
  body_ctrl_ = new BodyPostureCtrl<T>(robot, this);

  Test<T>::_state_list.push_back(body_up_ctrl_);
  Test<T>::_state_list.push_back(body_ctrl_);

  _sp = StateProvider<T>::getStateProvider();

  _SettingParameter();
  _folder_name = "/user/WBC_Controller/sim_data/";
  create_folder(_folder_name);
  printf("[Body Position Control Test] Constructed\n");
}

template <typename T>
BodyCtrlTest<T>::~BodyCtrlTest() {
  for (size_t i(0); i < Test<T>::_state_list.size(); ++i) {
    delete Test<T>::_state_list[i];
  }
}

template <typename T>
void BodyCtrlTest<T>::_TestInitialization() {
  // Yaml file name
  body_up_ctrl_->CtrlInitialization("CTRL_move_to_target_height");
  body_ctrl_->CtrlInitialization("CTRL_stance");
}

template <typename T>
int BodyCtrlTest<T>::_NextPhase(const int& phase) {
  int next_phase = phase + 1;
  printf("next phase: %i\n", next_phase);
  if (next_phase == NUM_BDCTRL_PHASE) {
    return BodyCtrlPhase::BDCTRL_body_ctrl;
  } else
    return next_phase;
}

template <typename T>
void BodyCtrlTest<T>::_SettingParameter() {
  typename std::vector<Controller<T>*>::iterator iter =
      Test<T>::_state_list.begin();

  while (iter < Test<T>::_state_list.end()) {
    if (Test<T>::_robot_type == RobotType::CHEETAH_3) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_body_ctrl_cheetah3.yaml");
      //
      ParamHandler handler(CheetahConfigPath "TEST_body_ctrl_cheetah3.yaml");
      handler.getBoolean("save_file", Test<T>::_b_save_file);
    } else if (Test<T>::_robot_type == RobotType::MINI_CHEETAH) {
      (*iter)->SetTestParameter(CheetahConfigPath
                                "TEST_body_ctrl_mini_cheetah.yaml");
      //
      ParamHandler handler(CheetahConfigPath
                           "TEST_body_ctrl_mini_cheetah.yaml");
    } else {
      printf("[Body Ctrl Test] Invalid robot type\n");
    }
    ++iter;
  }
}

template <typename T>
void BodyCtrlTest<T>::_UpdateTestOneStep() {
  // Data Save
  if (Test<T>::_b_save_file) {
    static int count(0);
    if (count % 10 == 0) {
      Vec3<T> body_ori =
          ori::quatToRPY(Test<T>::_robot->_state.bodyOrientation);
      saveVector(((BodyPostureCtrl<T>*)body_ctrl_)->_target_ori_command,
                 _folder_name, "cmd_body_ori_rpy");
      saveValue(_sp->_curr_time, _folder_name, "time");
      saveVector(body_ori, _folder_name, "body_ori_rpy");
      saveVector(_sp->_Q, _folder_name, "config");
      saveVector(_sp->_Qdot, _folder_name, "qdot");

      saveValue(Test<T>::_phase, _folder_name, "phase");
    }
    DVec<T> idx_jpos(cheetah::num_act_joint + 2);
    idx_jpos[0] = count;
    idx_jpos[1] = _sp->_curr_time * 1000.;
    idx_jpos.tail(cheetah::num_act_joint) =
        _sp->_Q.segment(6, cheetah::num_act_joint);
    saveVector(idx_jpos, _folder_name, "joint_pos");
    ++count;
  }
}

template <typename T>
void BodyCtrlTest<T>::_UpdateExtraData(Cheetah_Extra_Data<T>* ext_data) {
  ext_data->num_step = 0;
  ext_data->num_path_pt = 0;
}

template class BodyCtrlTest<double>;
template class BodyCtrlTest<float>;
