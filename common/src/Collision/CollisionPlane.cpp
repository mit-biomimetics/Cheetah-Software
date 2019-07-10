#include "Collision/CollisionPlane.h"

/*!
 * check whether the contact happens or not
 * cp_frame let you know which direction is normal (z) and which directions are
 * x and y w.r.t global frame. In the case of plane collition, the cp_frame is
 * always an identity matrix.
 */
template <typename T>
bool CollisionPlane<T>::ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                         Mat3<T>& cp_frame) {
  if (cp_pos[2] < _height) {
    penetration = cp_pos[2] - _height;
    cp_frame.setIdentity();
    return true;
  } else {
    return false;
  }
}

template class CollisionPlane<double>;
template class CollisionPlane<float>;
