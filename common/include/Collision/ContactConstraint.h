/*! @file ContactConstraint.h
 *  @brief Virtual class of Contact Constraint logic
 *
 * To make a child class, you need to implement virtual function,
 * UpdateExternalForces, UpdateQdot
 */

#ifndef CONTACT_CONSTRAINT_H
#define CONTACT_CONSTRAINT_H

#include <iostream>
#include "Collision/Collision.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Utilities/Utilities_print.h"
#include "cppTypes.h"

#define CC ContactConstraint<T>

template <typename T>
class ContactConstraint {
 public:
  ContactConstraint(FloatingBaseModel<T>* model)
      : _nContact(0), _nCollision(0) {
    _model = model;
    for (size_t i(0); i < _model->_nGroundContact; ++i) {
      _cp_force_list.push_back(Vec3<T>::Zero());
    }
  }

  virtual ~ContactConstraint() {}

  /*!
   * Add collision object
   * @param collision : collision objects
   */
  void AddCollision(Collision<T>* collision) {
    _collision_list.push_back(collision);
    ++_nCollision;
  }

  /*!
   * Used for spring damper based contact constraint method
   * @param K : Spring constant
   * @param D : Damping constant
   * @param dt : time step (sec)
   */
  virtual void UpdateExternalForces(T K, T D, T dt) = 0;

  /*!
   * Used for impulse based contact constraint method
   * @param state : full state of a floating system full state
   * @return state : update velocity
   */
  virtual void UpdateQdot(FBModelState<T>& state) = 0;

  /*!
   * For visualization
   * @return cp_pos_list : all contact point position in the global frame.
   */
  const vectorAligned<Vec3<T>>& getContactPosList() { return _cp_pos_list; }

  /*!
   * For visualization
   * @return cp_force_list : all linear contact force described in the global
   * frame.
   */
  const Vec3<T>& getGCForce(size_t idx) { return _cp_force_list[idx]; }

 protected:
  vectorAligned<Vec2<T>> deflectionRate;
  void _groundContactWithOffset(T K, T D);

  size_t _CheckContact();

  vectorAligned<Vec2<T>> _tangentialDeflections;

  size_t _nContact;
  size_t _nCollision;

  FloatingBaseModel<T>* _model;

  std::vector<Collision<T>*> _collision_list;
  std::vector<size_t> _idx_list;
  std::vector<T> _cp_resti_list;
  std::vector<T> _cp_mu_list;

  std::vector<T> _cp_penetration_list;
  vectorAligned<Vec3<T>> _cp_force_list;  // For All contact point w.r.t Global
  vectorAligned<Vec3<T>> _cp_local_force_list;
  vectorAligned<Vec3<T>> _cp_pos_list;
  vectorAligned<Mat3<T>> _cp_frame_list;
};

#endif
