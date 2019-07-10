/*! @file CollisionPlane.h
 *  @brief Collision logic for an infinite plane
 *
 *  Simplest collision, used for floor and global bounding box
 */

#ifndef COLLISIONPLANE_H
#define COLLISIONPLANE_H

#include <vector>
#include "Collision/Collision.h"
#include "cppTypes.h"

/*!
 * Class to represent infinite collision planes (like a flat ground).
 */
template <typename T>
class CollisionPlane : public Collision<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * Construct a new collision plane
   * @param location : coordinate transformation to collision plane (collision
   * surface is the xy-plane)
   * @param nContactPoints : number of contact points this collision plane will
   * need to handle.
   * @param mu : coefficient of friction
   * @param restitution  : rebounding ratio (v+/v-)
   * @param height  : height of this plane
   */
  CollisionPlane(const T& mu, const T& restitution, const T& height)
      : Collision<T>(mu, restitution), _height(height) {}

  virtual ~CollisionPlane() {}

  virtual bool ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                Mat3<T>& cp_frame);

 private:
  T _height;
};

#endif  // COLLISION_PLANE_H
