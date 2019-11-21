/*! @file Collision.h
 *  @brief Virtual class of Collision logic
 *
 * To make a child class, you need to implement virtual function,
 * ContactDetection, which checks if a point is in contact with the geometry
 */

#ifndef COLLISION_H
#define COLLISION_H

#include "cppTypes.h"

/*!
 * Abstract Collision Class
 */
template <typename T>
class Collision {
 public:
  /*!
   * Construct a new collision
   * @param mu : coefficient of friction
   * @param resti : coefficient of restitution (v_rebound / v_impact)
   */
  Collision(const T& mu, const T& resti) : _mu(mu), _restitution_coeff(resti) {}
  virtual ~Collision() {}

  /*!
   * virtual function for contact detection
   * @param cp_pos : contact point in the global frame
   * @param penetration : Size of the penetration to normal direction to the
   * collision object
   * @param cp_frame : Local frame that has normal axis (z) perpendicular the
   * contact surface
   */
  virtual bool ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                Mat3<T>& cp_frame) = 0;

  const T& getFrictionCoeff() { return _mu; }
  const T& getRestitutionCoeff() { return _restitution_coeff; }

 protected:
  T _mu, _restitution_coeff;
};

#endif  // COLLISION_H
