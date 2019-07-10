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

template <typename T>
class ActuatorModel {
 public:
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

  // compute actual actuator torque, given desired torque and speed.
  // takes into account friction (dry and damping), voltage limits, and torque
  // limits
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

  void setFriction(bool enabled) { _frictionEnabled = enabled; }

 private:
  T _gr, _kt, _R, _V, _damping, _dryFriction, _tauMax;
  bool _frictionEnabled = true;
};

#endif  // PROJECT_ACTUATORMODEL_H
