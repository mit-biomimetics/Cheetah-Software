/*! @file ActuatorModel.h
 *  @brief Model of actuator
 *  Includes friction, max torque, and motor torque speed curve.
 *
 *  The getTorque is used for torque at the joint, not torque at the motor.
 *  The provided frictions are for torques at the joint, not torque at the motor
 *  The R/KT are for the motor
 */

#ifndef PROJECT_ACTUATORMODEL_H
#define PROJECT_ACTUATORMODEL_H

#include "Utilities/utilities.h"

/*!
 * A model of an actuator containing friction and electrical effects
 */
template <typename T>
class ActuatorModel {
 public:

  /*!
   * Construct a new actuator model with the given parameters
   * @param gearRatio : Gear reduction
   * @param motorKT : Value of KT (torque constant) for the motor
   * @param motorR : Motor resistance
   * @param batteryV : Battery voltage
   * @param damping : Actuator damping (at the joint, Nm/(rad/sec))
   * @param dryFriction : Actuator dry friction (at the joint, Nm)
   * @param tauMax : Maximum torque output of the actuator
   */
  ActuatorModel(T gearRatio, T motorKT, T motorR, T batteryV, T damping,
                T dryFriction, T tauMax)
      : _gr(gearRatio),
        _kt(motorKT),
        _R(motorR),
        _V(batteryV),
        _damping(damping),
        _dryFriction(dryFriction),
        _tauMax(tauMax) {}

  ActuatorModel() {}

  // compute

  /*!
   * Compute actual actuator torque, given desired torque and speed.
   * takes into account friction (dry and damping), voltage limits, and torque
   * limits
   * @param tauDes : desired torque
   * @param qd : current actuator velocity (at the joint)
   * @return actual produced torque
   */
  T getTorque(T tauDes, T qd) {
    // compute motor torque
    T tauDesMotor = tauDes / _gr;        // motor torque
    T iDes = tauDesMotor / (_kt * 1.5);  // i = tau / KT
    // T bemf =  qd * _gr * _kt * 1.732;     // back emf
    T bemf = qd * _gr * _kt * 2.;       // back emf
    T vDes = iDes * _R + bemf;          // v = I*R + emf
    T vActual = coerce(vDes, -_V, _V);  // limit to battery voltage
    T tauActMotor =
        1.5 * _kt * (vActual - bemf) / _R;  // tau = Kt * I = Kt * V / R
    T tauAct = _gr * coerce(tauActMotor, -_tauMax, _tauMax);

    // add damping and dry friction
    if (_frictionEnabled)
      tauAct = tauAct - _damping * qd - _dryFriction * sgn(qd);

    return tauAct;
  }

  /*!
   * Control friction effects
   * @param enabled : enable/disable both dry and damping friction terms
   */
  void setFriction(bool enabled) { _frictionEnabled = enabled; }

 private:
  T _gr, _kt, _R, _V, _damping, _dryFriction, _tauMax;
  bool _frictionEnabled = true;
};

#endif  // PROJECT_ACTUATORMODEL_H
