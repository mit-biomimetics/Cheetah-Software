/*! @file DynamicsSimulator.cpp
 *  @brief Rigid Body Dynamics Simulator with Collisions
 *
 *  Combines ABA, Collisions, integrator, and any other external forces to run a
 * simulation. Doesn't do any graphics.
 */

#include "Dynamics/DynamicsSimulator.h"
#include "Collision/ContactImpulse.h"
#include "Collision/ContactSpringDamper.h"
#include "Utilities/Utilities_print.h"

/*!
 * Initialize the dynamics simulator by allocating memory for ABA matrices
 */
template <typename T>
DynamicsSimulator<T>::DynamicsSimulator(FloatingBaseModel<T> &model,
                                        bool useSpringDamper)
    : _model(model), _useSpringDamper(useSpringDamper) {
  if (_useSpringDamper) {
    _contact_constr = new ContactSpringDamper<T>(&_model);
  } else {
    _contact_constr = new ContactImpulse<T>(&_model);
  }

  _state.bodyVelocity = SVec<T>::Zero();
  _state.bodyPosition = Vec3<T>::Zero();
  _state.bodyOrientation = Quat<T>::Zero();
  _state.q = DVec<T>::Zero(_model._nDof - 6);
  _state.qd = DVec<T>::Zero(_model._nDof - 6);
  _dstate.qdd = DVec<T>::Zero(_model._nDof - 6);
  _lastBodyVelocity.setZero();
}

/*!
 * Take one simulation step
 * @param dt : timestep duration
 * @param tau : joint torques
 */
template <typename T>
void DynamicsSimulator<T>::step(T dt, const DVec<T> &tau, T kp, T kd ) {
  // fwd-kin on gc points
  // compute ground contact forces
  // aba
  // integrate
  forwardKinematics();           // compute forward kinematics
  updateCollisions(dt, kp, kd);  // process collisions
  // Process Homing
  if( _homing.active_flag) {

    Mat3<T> R10_des = rpyToRotMat(_homing.rpy);              // R10_des
    Mat3<T> R10_act = _model.getOrientation(5).transpose();  // R10
    Mat3<T> eR01 = R10_des.transpose()*R10_act;              // eR * R01 = R01_des
    
    Vec4<T> equat = rotationMatrixToQuaternion(eR01.transpose());
    Vec3<T> angle_axis = quatToso3(equat); // in world frame

    Vec3<T> p = _model.getPosition(5);
    Vec3<T> f = _homing.kp_lin*(_homing.position - p)-_homing.kd_lin*_model.getLinearVelocity(5);

    // Note: External forces are spatial forces in the {0} frame. 
    _model._externalForces.at(5) += forceToSpatialForce(f,p);
    _model._externalForces.at(5).head(3) += _homing.kp_ang*angle_axis - _homing.kd_ang*_model.getAngularVelocity(5);

  }


  runABA(tau);                   // dynamics algorithm
  integrate(dt);                 // step forward

  _model.setState(_state);
  _model.resetExternalForces();  // clear external forces
  _model.resetCalculationFlags();
}

/*!
 * Run spring-damper collisions
 * @param dt : timestep
 * @param kp : spring constant
 * @param kd : damping constant
 */
template <typename T>
void DynamicsSimulator<T>::updateCollisions(T dt, T kp, T kd) {
  _model.forwardKinematics();
  _contact_constr->UpdateExternalForces(kp, kd, dt);
}

/*!
 * Integrate the floating base state
 * @param dt timestep
 */
template <typename T>
void DynamicsSimulator<T>::integrate(T dt) {
  if (_useSpringDamper) {
    Vec3<T> omegaBody = _state.bodyVelocity.template block<3, 1>(0, 0);
    Mat6<T> X = createSXform(quaternionToRotationMatrix(_state.bodyOrientation),
                             _state.bodyPosition);
    RotMat<T> R = rotationFromSXform(X);
    Vec3<T> omega0 = R.transpose() * omegaBody;

    // actual integration
    _state.qd += _dstate.qdd * dt;
    _state.q += _state.qd * dt;

    _state.bodyVelocity += _dstate.dBodyVelocity * dt;
    _state.bodyPosition += _dstate.dBodyPosition * dt;
    _state.bodyOrientation = integrateQuat(_state.bodyOrientation, omega0, dt);
  } else {
    // actual integration
    // Velocity Update by integrating acceleration
    _state.qd += _dstate.qdd * dt;
    _state.bodyVelocity += _dstate.dBodyVelocity * dt;

    // Contact Constraint Velocity Updated
    _contact_constr->UpdateQdot(_state);

    // Prepare body velocity integration
    RotMat<T> R_body = quaternionToRotationMatrix(_state.bodyOrientation);

    _dstate.dBodyPosition =
        R_body.transpose() * _state.bodyVelocity.template block<3, 1>(3, 0);
    Vec3<T> omegaBody = _state.bodyVelocity.template block<3, 1>(0, 0);

    // Position Update
    _state.q += _state.qd * dt;
    _state.bodyPosition += _dstate.dBodyPosition * dt;
    _state.bodyOrientation =
        integrateQuatImplicit(_state.bodyOrientation, omegaBody, dt);
    _dstate.dBodyVelocity = (_state.bodyVelocity - _lastBodyVelocity) / dt;
    _lastBodyVelocity = _state.bodyVelocity;
  }
}

template class DynamicsSimulator<double>;

template class DynamicsSimulator<float>;
