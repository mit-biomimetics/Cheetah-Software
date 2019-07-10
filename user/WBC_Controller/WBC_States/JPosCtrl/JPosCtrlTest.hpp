#ifndef JPOS_CTR_TEST_H
#define JPOS_CTR_TEST_H

#include <WBC_States/Test.hpp>

enum JPosCtrlPhase {
  JPCTRL_move_to_target = 0,
  JPCTRL_stay = 1,
  JPCTRL_swing = 2,
  NUM_JPCTRL_PHASE
};

template <typename T>
class JPosCtrlTest : public Test<T> {
 public:
  JPosCtrlTest(FloatingBaseModel<T>*, const RobotType&);
  virtual ~JPosCtrlTest();

 protected:
  virtual void _TestInitialization();
  virtual int _NextPhase(const int& phase);
  virtual void _UpdateExtraData(Cheetah_Extra_Data<T>* ext_data);
  void _SettingParameter();

  Controller<T>* _ini_jpos_ctrl;
  Controller<T>* _jpos_stay;
  Controller<T>* _jpos_swing;
};

#endif
