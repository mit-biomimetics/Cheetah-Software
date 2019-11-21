#ifndef WBC_CONTROLLER_H
#define WBC_CONTROLLER_H

#include <FSM_States/ControlFSMData.h>
#include <Dynamics/FloatingBaseModel.h>
#include <Dynamics/Quadruped.h>
#include "cppTypes.h"
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBIC/KinWBC.hpp>

#include <lcm-cpp.hpp>
#include "wbc_test_data_t.hpp"

#define WBCtrl WBC_Ctrl<T>

class MIT_UserParameters;

template<typename T>
class WBC_Ctrl{
  public:
    WBC_Ctrl(FloatingBaseModel<T> model);
    virtual ~WBC_Ctrl();

    void run(void * input, ControlFSMData<T> & data);
    void setFloatingBaseWeight(const T & weight){
      _wbic_data->_W_floating = DVec<T>::Constant(6, weight);
    }

  protected:
    virtual void _ContactTaskUpdate(void * input, ControlFSMData<T> & data) = 0;
    virtual void _ContactTaskUpdateTEST(void * input, ControlFSMData<T> & data){
      (void)input;
      (void)data;
    }
    virtual void _LCM_PublishData(){}
    void _UpdateModel(const StateEstimate<T> & state_est, const LegControllerData<T> * leg_data);
    void _UpdateLegCMD(ControlFSMData<T> & data);
    void _ComputeWBC();

    KinWBC<T>* _kin_wbc;
    WBIC<T>* _wbic;
    WBIC_ExtraData<T>* _wbic_data;

    FloatingBaseModel<T> _model;
    std::vector<ContactSpec<T> * > _contact_list;
    std::vector<Task<T> * > _task_list;

    DMat<T> _A;
    DMat<T> _Ainv;
    DVec<T> _grav;
    DVec<T> _coriolis;

    FBModelState<T> _state;

    DVec<T> _full_config;
    DVec<T> _tau_ff;
    DVec<T> _des_jpos;
    DVec<T> _des_jvel;

    std::vector<T> _Kp_joint, _Kd_joint;
    //std::vector<T> _Kp_joint_swing, _Kd_joint_swing;

    unsigned long long _iter;

    lcm::LCM _wbcLCM;
    wbc_test_data_t _wbc_data_lcm;
};
#endif
