/*! @file CollisionMesh.h
 *  @brief Collision logic for a mesh
 */

#ifndef COLLISION_MESH_H
#define COLLISION_MESH_H

#include <vector>
#include "Collision/Collision.h"
#include "Utilities/utilities.h"
#include "cppTypes.h"

/*!
 * Class to represent box collision
 */
template <typename T>
class CollisionMesh : public Collision<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * Construct a new collision Mesh
   * @param mu : coefficient of friction
   * @param restitution  : rebounding ratio (v+/v-)
   * @param grid : grid size (meter) of the given height map
   * @param left_corner_loc : global location of left bottom edge of the height
   * map
   * @param height_map : Height map of mesh
   */
  CollisionMesh(const T& mu, const T& restitution, const T& grid,
                const Vec3<T>& left_corner_loc, const DMat<T>& height_map)
      : Collision<T>(mu, restitution),
        _left_corner_loc(left_corner_loc),
        _height_map(height_map),
        _grid(grid) {
    _x_max = (_height_map.rows() - 1) * _grid;
    _y_max = (_height_map.cols() - 1) * _grid;
  }

  virtual ~CollisionMesh() {}
  virtual bool ContactDetection(const Vec3<T>& cp_pos, T& penetration,
                                Mat3<T>& cp_frame);

 private:
  T _size[3];

  Vec3<T> _left_corner_loc;
  DMat<T> _height_map;

  T _grid;
  T _x_max;
  T _y_max;
};

#endif  // COLLISION_MESH_H
