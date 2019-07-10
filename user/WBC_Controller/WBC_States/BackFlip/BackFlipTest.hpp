#ifndef BACK_FLIP_TEST_Cheetah
#define BACK_FLIP_TEST_Cheetah

#include <Dynamics/FloatingBaseModel.h>
#include <WBC_States/Test.hpp>
#include "DataReader.hpp"

template <typename T>
class StateProvider;

enum BFlipPhase {
  BFlip_body_up_ctrl = 0,
  // BFlip_backflip_pre = 1,
  BFlip_backflip = 1,  // 2,
  BFlip_backflip_landing = 2,
  BFlip_backflip_end = 3,  // 3,
  NUM_BFlip_PHASE = 4,
};

template <typename T>
class BackFlipTest : public Test<T> {
 public:
  BackFlipTest(FloatingBaseModel<T>*, const RobotType&);
  virtual ~BackFlipTest();

 protected:
  virtual void _UpdateTestOneStep();
  virtual void _TestInitialization();
  virtual void _UpdateExtraData(Cheetah_Extra_Data<T>* ext_data);
  virtual int _NextPhase(const int& phase);
  void _SettingParameter();

  std::string _folder_name;

  Controller<T>* body_up_ctrl_;
  // Controller<T>* backflip_pre_;
  Controller<T>* backflip_ctrl_;
  Controller<T>* backflip_landing_;
  Controller<T>* body_ctrl_;

  DataReader* _data_reader;
  StateProvider<T>* _sp;
};

#endif
