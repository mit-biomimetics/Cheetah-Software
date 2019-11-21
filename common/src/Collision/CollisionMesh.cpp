/*!
 * @file CollisionMesh.cpp
 * @brief Collision logic for a mesh/height map
 */

#include "Collision/CollisionMesh.h"
#include "Utilities/Utilities_print.h"

/*!
 * check whether the contact happens or not
 * cp_frame let you know which direction is normal (z) and which directions are
 * x and y. The frame is basically contact coordinate w.r.t global.
 */
template <typename T>
bool CollisionMesh<T>::ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                        Mat3<T>& cp_frame) {
  // Contact detection
  Vec3<T> cp_pos_in_height_map = cp_pos - _left_corner_loc;
  if ((0 < cp_pos_in_height_map[0]) && (cp_pos_in_height_map[0] < _x_max)) {
    if ((0 < cp_pos_in_height_map[1]) && (cp_pos_in_height_map[1] < _y_max)) {
      int x_idx = floor(cp_pos_in_height_map[0] / _grid);
      int y_idx = floor(cp_pos_in_height_map[1] / _grid);

      Vec3<T> vec1;
      vec1[0] = _grid;
      vec1[1] = _grid;
      vec1[2] = _height_map(x_idx + 1, y_idx + 1) - _height_map(x_idx, y_idx);

      Vec3<T> vec2;
      vec2[0] = _grid;
      vec2[1] = -_grid;
      vec2[2] = _height_map(x_idx + 1, y_idx) - _height_map(x_idx, y_idx + 1);

      vec1.normalize();
      vec2.normalize();

      Vec3<T> normal = vec2.cross(vec1);
      // vec2 is x axis
      Vec3<T> y_axis = normal.cross(vec2);

      cp_frame.template block<3, 1>(0, 0) = vec2;
      cp_frame.template block<3, 1>(0, 1) = y_axis;
      cp_frame.template block<3, 1>(0, 2) = normal;

      // Vector pointing middle of four points
      T height_ave = 0.25 * _height_map(x_idx, y_idx) +
                     0.25 * _height_map(x_idx, y_idx + 1) +
                     0.25 * _height_map(x_idx + 1, y_idx) +
                     0.25 * _height_map(x_idx + 1, y_idx + 1);

      Vec3<T> middle_pt;
      middle_pt[0] = _grid * x_idx + 0.5 * _grid;
      middle_pt[1] = _grid * y_idx + 0.5 * _grid;
      middle_pt[2] = height_ave;

      Vec3<T> local_diff =
          cp_frame.transpose() * (cp_pos_in_height_map - middle_pt);

      if (local_diff[2] < 0.) {
        penetration = local_diff[2];
        return true;
      }
    }
  }

  return false;
}

template class CollisionMesh<double>;
template class CollisionMesh<float>;
