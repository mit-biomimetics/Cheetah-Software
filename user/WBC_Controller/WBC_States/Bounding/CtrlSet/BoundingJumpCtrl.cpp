#include "BoundingJumpCtrl.hpp"
#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <WBC_States/Bounding/BoundingTest.hpp>
#include <Math/orientation_tools.h>
#include <Utilities/Utilities_print.h>

#include <WBC/WBIC/WBIC.hpp>
#include <WBC/WBLC/KinWBC.hpp>
#include <Utilities/Timer.h>

#include <WBC_States/Cheetah_DynaCtrl_Definition.h>
#include <WBC_States/Bounding/TaskSet/BodyRyRzTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalHeadPosTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalPosTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalRollTask.hpp>
#include <WBC_States/Bounding/TaskSet/LocalTailPosTask.hpp>
#include <WBC_States/StateProvider.hpp>
#include <WBC_States/common/ContactSet/SingleContact.hpp>


template <typename T>
BoundingJumpCtrl<T>::BoundingJumpCtrl(BoundingTest<T>* bounding_test,
                                    const FloatingBaseModel<T>* robot)
  : Controller<T>(robot),
  _ctrl_start_time(0.),
  _des_jpos(cheetah::num_act_joint),
  _des_jvel(cheetah::num_act_joint),
  _des_jacc(cheetah::num_act_joint),
  _bounding_test(bounding_test){

  _local_roll_task = new LocalRollTask<T>(Ctrl::_robot_sys);
  _body_ryrz_task = new BodyRyRzTask<T>(Ctrl::_robot_sys);

  _fr_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::FR, linkID::FR_abd);
  _fl_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::FL, linkID::FL_abd);
  _hr_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::HR, linkID::HR_abd);
  _hl_foot_local_task =
      new LocalPosTask<T>(Ctrl::_robot_sys, linkID::HL, linkID::HL_abd);

  _local_head_pos_task = new LocalHeadPosTask<T>(Ctrl::_robot_sys);
  _local_tail_pos_task = new LocalTailPosTask<T>(Ctrl::_robot_sys);

  _fr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FR);
  _fl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::FL);
  _hr_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HR);
  _hl_contact = new SingleContact<T>(Ctrl::_robot_sys, linkID::HL);

  _kin_wbc = new KinWBC<T>(cheetah::dim_config);
  _wbic = new WBIC<T>(cheetah::dim_config, &(Ctrl::_contact_list), &(Ctrl::_task_list));

  _wbic_data = new WBIC_ExtraData<T>();
  _wbic_data->_W_floating = DVec<T>::Constant(6, 5.);
  _wbic_data->_W_rf = DVec<T>::Constant(6, 1.0);
  _wbic_data->_W_rf[2] = 50.0;
  _wbic_data->_W_rf[5] = 50.0;

  _sp = StateProvider<T>::getStateProvider();
}

template<typename T>
BoundingJumpCtrl<T>::~BoundingJumpCtrl(){
  delete _fr_foot_local_task;
  delete _fl_foot_local_task;
  delete _hr_foot_local_task;
  delete _hl_foot_local_task;

  delete _local_head_pos_task;
  delete _local_tail_pos_task;

  delete _fr_contact;
  delete _fl_contact; 
  delete _hr_contact; 
  delete _hl_contact; 
}

template<typename T>
void BoundingJumpCtrl<T>::_optimize_jump(const DVec<T> & X0, T xdot_target, T zdot_target ){

  int N = 31;  // number of control intervals
  int N_fc = 15;
  auto opti = casadi::Opti();  // Optimization problem

  Slice all;
  // ---- decision variables ---------
  auto X = opti.variable(6, N+1);  // state trajectory
  auto P = opti.variable(4);  // Landing Location
  auto F = opti.variable(2, N);  // Reaction force (15: front, 1: air, 15: hind)
  auto cost = opti.variable();

  //DVec<double> X0(6); X0.setZero();
  //X0<< 0.0, 0.2, 0.2, 0.5, -0.1, 0.1;

  DVec<double> Xf(6); Xf.setZero();
  //Xf<< 0.1, 0.2, -M_PI/4., 0.1, 0.3, 0.0;
  Xf<< 0.5, 0.23, -M_PI/4., xdot_target, zdot_target, 0.1;

  cost = 
    (X(0,N) - Xf(0)) * (X(0,N) - Xf(0))
    + (X(1,N) - Xf(1)) * (X(1,N) - Xf(1))
    + (X(2,N) - Xf(2)) * (X(2,N) - Xf(2))
    + (X(3,N) - Xf(3)) * (X(3,N) - Xf(3))
    + (X(4,N) - Xf(4)) * (X(4,N) - Xf(4))
    + (X(5,N) - Xf(5)) * (X(5,N) - Xf(5));
  //
  // ---- objective          ---------
  opti.minimize(cost);  // race in minimal time

  // ---- dynamic constraints --------
  auto dt = 0.01;
  double mu(0.5);

  for (int k = 0; k < N; ++k) {
    auto dX = simple_dyn(X(all, k), F(all, k), P, k);
    auto X_next = X(all, k) + dt * dX;
    opti.subject_to(X(all, k + 1) == X_next);  // close the gaps

    // ---- contact constraints -----------
    opti.subject_to(-mu * F(1,k) <= F(0, k) <= mu* F(1,k)); 
    opti.subject_to(0<=F(1,k) <= 300.); 
    // Kinematics constraints
    opti.subject_to (-0.15 <= P(0) <= 0.15);
    opti.subject_to (-0.3 <= P(2) <= -0.1);
    opti.subject_to (P(1) == 0.);
    opti.subject_to (P(3) == 0.);
  }
  opti.subject_to(F(1, N_fc) == 0.);

  for(int i(0); i<6; ++i){
    opti.subject_to(X(i, 1) == X0(i)); 
  }
  // ---- solve NLP   
  Dict plugin_opts;
  plugin_opts["ipopt.print_level"] = 0;
  Dict solver_opts;
  opti.solver("ipopt", plugin_opts, solver_opts);     // set numerical backend
  auto sol = opti.solve();  // actual solve

  _xpos_list = std::vector<T>(sol.value(X(0, all) ) );
  _zpos_list = std::vector<T>(sol.value(X(1, all) ) );
  _theta_list = std::vector<T>(sol.value(X(2, all) ) );

  _xvel_list = std::vector<T>(sol.value(X(3, all) ) );
  _zvel_list = std::vector<T>(sol.value(X(4, all) ) );
  _theta_dot_list = std::vector<T>(sol.value(X(5, all) ) );

  _fr_x_list = std::vector<T>(sol.value(F(0, all)));
  _fr_z_list = std::vector<T>(sol.value(F(1, all)));

  _P_foot = std::vector<T>(sol.value(P));
}

template<typename T>
MX BoundingJumpCtrl<T>::simple_dyn(const MX& x, const MX& u, const MX& P, int i) { 
  T I(0.21);
  T M(8.91);
  T g(9.81);

  if( i<15 ){
    return vertcat (vertcat(x(3), x(4), x(5)), 
        vertcat(u(0)/M, u(1)/M -g), ((P(1) - x(1))*u(0) - (P(0) - x(0)) * u(1))/I );  // Velocity
  }else{
    return vertcat (vertcat(x(3), x(4), x(5)), 
        vertcat(u(0)/M, u(1)/M -g), ((P(3) - x(1))*u(0) - (P(2) - x(0)) * u(1))/I );  // Velocity
  }
}

template <typename T>
void BoundingJumpCtrl<T>::OneStep(void* _cmd) {
  Ctrl::_PreProcessing_Command();

  if(_iter%10 == 0) ++_motion_count;
  
  // Update Time
  Ctrl::_state_machine_time = _sp->_curr_time - _ctrl_start_time;
  DVec<T> gamma = DVec<T>::Zero(cheetah::num_act_joint);

  _contact_update();
  _update_motion_cmd();
  _compute_torque_wbic(gamma);
  // Update Contact
  for (size_t leg(0); leg < cheetah::num_leg; ++leg) {
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
      ((LegControllerCommand<T>*)_cmd)[leg].tauFeedForward[jidx] =
          gamma[cheetah::num_leg_joint * leg + jidx];

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
  ++_iter;

  Ctrl::_PostProcessing_Command();
}

template <typename T>
void BoundingJumpCtrl<T>::_compute_torque_wbic(DVec<T>& gamma) {
  _kin_wbc->FindConfiguration(_sp->_Q, Ctrl::_task_list, Ctrl::_contact_list,
                              _des_jpos, _des_jvel, _des_jacc);

  // WBIC
  _wbic->UpdateSetting(Ctrl::_A, Ctrl::_Ainv, Ctrl::_coriolis, Ctrl::_grav);
  _wbic->MakeTorque(gamma, _wbic_data);
}


template<typename T>
void BoundingJumpCtrl<T>::_contact_update(){
  Ctrl::_contact_list.clear();

  if(_motion_count < 15){
    Ctrl::_contact_list.push_back(_fr_contact);
    Ctrl::_contact_list.push_back(_fl_contact);
    _fr_contact->UpdateContactSpec();
    _fl_contact->UpdateContactSpec();
  }else{
    Ctrl::_contact_list.push_back(_hr_contact);
    Ctrl::_contact_list.push_back(_hl_contact);
    _hr_contact->UpdateContactSpec();
    _hl_contact->UpdateContactSpec();
  }
}

template<typename T>
void BoundingJumpCtrl<T>::_update_motion_cmd(){
  Ctrl::_task_list.clear();
  Ctrl::_task_list.push_back(_local_roll_task);
  Ctrl::_task_list.push_back(_body_ryrz_task);

  if(_motion_count < 15){ // Front foot contact
    Ctrl::_task_list.push_back(_local_head_pos_task);
  }else{

  }

}
template <typename T>
void BoundingJumpCtrl<T>::SetTestParameter(const std::string &test_file) {
  _param_handler = new ParamHandler(test_file);
  
  // Joint level feedback gain
  _param_handler->getVector<T>("Kp_joint", _Kp_joint);
  _param_handler->getVector<T>("Kd_joint", _Kd_joint);
}

template<typename T>
void BoundingJumpCtrl<T>::CtrlInitialization(const std::string &category_name) {
  _param_handler->getValue<T>(category_name, "default_z_vel", _z_vel_target);
}


template <typename T>
bool BoundingJumpCtrl<T>::EndOfPhase() {
  if(_iter > 298){
    return true;
  }
  return false;
}

template <typename T>
void BoundingJumpCtrl<T>::FirstVisit() {
  _iter = 0;
  _motion_count = -1;

  Vec3<T> front_foot = 0.5*Ctrl::_robot_sys->_pGC[linkID::FR] + 0.5*Ctrl::_robot_sys->_pGC[linkID::FL];
  Vec3<T> body_pos = Ctrl::_robot_sys->_state.bodyPosition;
  Vec3<T> ori_rpy = ori::quatToRPY(Ctrl::_robot_sys->_state.bodyOrientation);
  Vec3<T> local_pos = body_pos - front_foot;

  Vec3<T> front_foot_vel = 0.5*Ctrl::_robot_sys->_vGC[linkID::FR] + 0.5*Ctrl::_robot_sys->_vGC[linkID::FL];
  SVec<T> body_vel = Ctrl::_robot_sys->_state.bodyVelocity;
  Vec3<T> local_body_vel = body_vel.tail(3) - front_foot_vel;
  
  pretty_print(front_foot_vel, std::cout, "front foot vel");
  pretty_print(body_pos, std::cout, "body pos");
  pretty_print(ori_rpy, std::cout, "ori rpy");

  DVec<T> X0(6);
  X0 << local_pos[0], local_pos[2], ori_rpy[1], local_body_vel[0], local_body_vel[2], body_vel[1];
  pretty_print(X0, std::cout, "X0");
  Timer timer;
  _optimize_jump(X0, 0.5, 0.5);
  printf("computation time: %f\n", timer.getMs());
}


template <typename T>
void BoundingJumpCtrl<T>::LastVisit() {}

 
template class BoundingJumpCtrl<double>;
template class BoundingJumpCtrl<float>;
