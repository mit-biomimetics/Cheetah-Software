#include "JPosCtrlTest.hpp"

#include <ParamHandler/ParamHandler.hpp>
#include <WBC_States/common/CtrlSet/JPosCtrl.hpp>

template <typename T>
JPosCtrlTest<T>::JPosCtrlTest(FloatingBaseModel<T>* robot,
                              const RobotType& robot_type)
    : Test<T>(robot, robot_type) {
  Test<T>::_phase = JPosCtrlPhase::JPCTRL_move_to_target;
  Test<T>::_state_list.clear();

  _ini_jpos_ctrl = new JPosCtrl<T>(robot);
  _jpos_stay = new JPosCtrl<T>(robot);
  _jpos_swing = new JPosCtrl<T>(robot);

  Test<T>::_state_list.push_back(_ini_jpos_ctrl);
  Test<T>::_state_list.push_back(_jpos_stay);
  Test<T>::_state_list.push_back(_jpos_swing);

  _SettingParameter();

  printf("[Joint Position Control Test] Constructed\n");
}

template <typename T>
JPosCtrlTest<T>::~JPosCtrlTest() {
  for (size_t i(0); i < Test<T>::_state_list.size(); ++i) {
    delete Test<T>::_state_list[i];
  }
}

template <typename T>
void JPosCtrlTest<T>::_TestInitialization() {
  // Yaml file name
  _ini_jpos_ctrl->CtrlInitialization("CTRL_jpos_move_to_target");
  _jpos_stay->CtrlInitialization("CTRL_jpos_stay");
  _jpos_swing->CtrlInitialization("CTRL_jpos_swing");
}

template <typename T>
int JPosCtrlTest<T>::_NextPhase(const int& phase) {
  int next_phase = phase + 1;
  printf("next phase: %i\n", next_phase);
  if (next_phase == NUM_JPCTRL_PHASE) {
    return JPosCtrlPhase::JPCTRL_stay;
  } else
    return next_phase;
}

template <typename T>
void JPosCtrlTest<T>::_SettingParameter() {
  typename std::vector<Controller<T>*>::iterator iter =
      Test<T>::_state_list.begin();
  while (iter < Test<T>::_state_list.end()) {
    (*iter)->SetTestParameter(CheetahConfigPath "TEST_jpos_ctrl.yaml");
    ++iter;
  }
}

template <typename T>
void JPosCtrlTest<T>::_UpdateExtraData(Cheetah_Extra_Data<T>* ext_data) {
  ext_data->num_step = 0;
  ext_data->num_path_pt = 0;
}

template class JPosCtrlTest<double>;
template class JPosCtrlTest<float>;
