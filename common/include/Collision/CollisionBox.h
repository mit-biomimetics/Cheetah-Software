/*! @file CollisionBox.h
 *  @brief Collision logic for a box
 */

#ifndef COLLISION_BOX_H
#define COLLISION_BOX_H

#include <vector>
#include "Collision/Collision.h"
#include "Utilities/utilities.h"
#include "cppTypes.h"

/*!
 * Class to represent box collision
 */
template <typename T>
class CollisionBox : public Collision<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * Construct a new collision box
   * @param mu : coefficient of friction
   * @param restitution  : rebounding ratio (v+/v-)
   * @param depth : depth of the box
   * @param width : width of the box
   * @param height : height of the box
   * @param position : position of the box w.r.t the global frame
   * @param ori : orientation of the box w.r.t the global frame
   */
  CollisionBox(const T& mu, const T& restitution, const T& depth,
               const T& width, const T& height, const Vec3<T>& position,
               const Mat3<T>& ori)
      : Collision<T>(mu, restitution), _position(position), _orientation(ori) {
    _size[0] = depth;
    _size[1] = width;
    _size[2] = height;
  }

  virtual ~CollisionBox() {}
  virtual bool ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                Mat3<T>& cp_frame);

 private:
  T _size[3];

  Vec3<T> _position;
  Mat3<T> _orientation;
};

#endif  // COLLISION_BOX_H
