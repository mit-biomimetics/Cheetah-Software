#ifndef LOCOMOTION_CONTROLLER
#define LOCOMOTION_CONTROLLER


#include <WBC_Ctrl/WBC_Ctrl.hpp>

template<typename T>
class LocomotionCtrlData{
  public:
    Vec3<T> pBody_des;
    Vec3<T> vBody_des;
    Vec3<T> aBody_des;
    Vec3<T> pBody_RPY_des;
    Vec3<T> vBody_Ori_des;

    Vec3<T> pFoot_des[4];
    Vec3<T> vFoot_des[4];
    Vec3<T> aFoot_des[4];
    Vec3<T> Fr_des[4];

    Vec4<T> contact_state;
};

template<typename T>
class LocomotionCtrl: public WBC_Ctrl<T>{
  public:
    LocomotionCtrl(FloatingBaseModel<T> model);
    virtual ~LocomotionCtrl();

  protected:
    virtual void _ContactTaskUpdate(
        void * input, ControlFSMData<T> & data);
    virtual void _ContactTaskUpdateTEST(void * input, ControlFSMData<T> & data);
    void _ParameterSetup(const MIT_UserParameters* param);
    void _CleanUp();
    virtual void _LCM_PublishData();

    LocomotionCtrlData<T>* _input_data;

    Task<T>* _body_pos_task;
    Task<T>* _body_ori_task;

    Task<T>* _foot_task[4];
    ContactSpec<T>* _foot_contact[4];

    Vec3<T> pre_foot_vel[4];

    Vec3<T> _Fr_result[4];
    Quat<T> _quat_des;
};

#endif

