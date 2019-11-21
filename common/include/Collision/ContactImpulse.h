/*!
 * @file ContactImpulse.h
 * @brief Implementation of impulse-based contact dynamics
 *
 * These are the default contact dynamics of the simulator.
 */

#ifndef CONTACT_IMPULSE_H
#define CONTACT_IMPULSE_H

#include "ContactConstraint.h"
#include "Dynamics/FloatingBaseModel.h"

/*!
 * Impulse-based contact dynamics for a FloatingBaseModel
 */
template <typename T>
class ContactImpulse : public ContactConstraint<T> {
 public:

  /*!
   * Construct a new contact dynamics model for a given FloatingBaseModel
   * @param model : The floating base model used for contacts
   */
  ContactImpulse(FloatingBaseModel<T>* model) : ContactConstraint<T>(model) {
    _iter_lim = 10;
    _nDof = CC::_model->_nDof;
    _b_debug = false;
  }

  virtual ~ContactImpulse() {}

  /*!
   * Impulse-based contact dynamics do not adjust the external forces, they modify q_dot
   * instead.
   *
   * Used for spring damper based contact constraint method,
   * therefore almost empty function in this class.
   * Only used to update dt and penetration recovery ratio
   * @param dt : time step (sec)
   */
  virtual void UpdateExternalForces(T K, T D, T dt) {
    (void)K;
    (void)D;
    // Set penetration recovery ration here
    _penetration_recover_ratio = 0.0 / dt;
    _dt = dt;
  }

  virtual void UpdateQdot(FBModelState<T>& state);

 protected:
  size_t _iter_lim;
  bool _b_debug;
  T _tol = 1.e-6;
  T _penetration_recover_ratio;
  size_t _nDof;
  void _UpdateVelocity(DVec<T>& qdot);
  void _UpdateQdotOneDirection(size_t idx,
                               const vectorAligned<D3Mat<T> >& Jc_list,
                               const T* lambda_list,
                               const vectorAligned<DVec<T> > AinvB_list,
                               const T* des_vel_list, const T* min_list,
                               const T* max_list, DVec<T>& qdot);
  size_t _iter_sum;

 private:
  T _dt;
};

#endif
