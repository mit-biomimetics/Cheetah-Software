/*! @file VisualizationData.h
 *  @brief Data sent from robot code to simulator GUI for debugging
 *
 */

#ifndef VISUALIZATION_DATA_H
#define VISUALIZATION_DATA_H

#define VISUALIZATION_MAX_PATH_POINTS 2000
#define VISUALIZATION_MAX_PATHS 10
#define VISUALIZATION_MAX_ITEMS 10000
#define VISUALIZATION_MAX_CHEETAHS 0

#define VISUALIZATION_MAX_MESHES 5
#define VISUALIZATION_MAX_MESH_GRID 150

#include "cppTypes.h"

/*!
 * Debugging sphere
 */
struct SphereVisualization {
  Vec3<float> position;
  Vec4<float> color;
  double radius;
};

/*!
 * Debugging box
 */
struct BlockVisualization {
  Vec3<float> dimension;
  Vec3<float> corner_position;
  Vec3<float> rpy;
  Vec4<float> color;
};

/*!
 * Debugging arrow
 */
struct ArrowVisualization {
  Vec3<float> base_position;
  Vec3<float> direction;
  Vec4<float> color;
  float head_width;
  float head_length;
  float shaft_width;
};

/*!
 * Debugging robot (draws the same type of robot as currently simulating)
 */
struct CheetahVisualization {
  Vec12<float> q;
  Quat<float> quat;
  Vec3<float> p;
  Vec4<float> color;
};

/*!
 * Debugging "path"
 */
struct PathVisualization {
  size_t num_points = 0;
  Vec4<float> color;
  Vec3<float> position[VISUALIZATION_MAX_PATH_POINTS];
  void clear() {
    num_points = 0;
  }
};

/*!
 * Debugging Cone
 */
struct ConeVisualization {
  Vec3<float> point_position;
  Vec3<float> direction;
  Vec4<float> color;
  double radius;
};

/*!
 * Mesh Visualization
 */
struct MeshVisualization {
  Vec3<float> left_corner;
  Eigen::Matrix<float, VISUALIZATION_MAX_MESH_GRID, VISUALIZATION_MAX_MESH_GRID> height_map;

  int rows, cols;

  float grid_size;
  float height_max;
  float height_min;
};


/*!
 * Collection of all debugging data
 */
struct VisualizationData {
  size_t num_paths = 0, num_arrows = 0, num_cones = 0, num_spheres = 0,
         num_blocks = 0, num_meshes = 0;
  SphereVisualization spheres[VISUALIZATION_MAX_ITEMS];
  BlockVisualization blocks[VISUALIZATION_MAX_ITEMS];
  ArrowVisualization arrows[VISUALIZATION_MAX_ITEMS];
  ConeVisualization cones[VISUALIZATION_MAX_ITEMS];
  PathVisualization paths[VISUALIZATION_MAX_PATHS];
  MeshVisualization meshes[VISUALIZATION_MAX_MESHES];

  /*!
   * Remove all debug data
   */
  void clear() {
    num_paths = 0, num_arrows = 0, num_cones = 0, num_spheres = 0, num_blocks = 0;
    num_meshes = 0;
  }

  /*!
   * Add a new sphere
   * @return A sphere, or nullptr if there isn't enough room
   */
  SphereVisualization* addSphere() {
    if(num_spheres < VISUALIZATION_MAX_ITEMS) {
      return &spheres[num_spheres++];
    }
    return nullptr;
  }

  /*!
   * Add a new box
   * @return A box, or nullptr if there isn't enough room
   */
  BlockVisualization* addBlock() {
    if(num_blocks < VISUALIZATION_MAX_ITEMS) {
      return &blocks[num_blocks++];
    }
    return nullptr;
  }

  /*!
   * Add a new arrow
   * @return An arrow, or nullptr if there isn't enough room
   */
  ArrowVisualization* addArrow() {
    if(num_arrows < VISUALIZATION_MAX_ITEMS) {
      return &arrows[num_arrows++];
    }
    return nullptr;
  }

  /*!
   * Add a new cone
   * @return A cone, or nullptr if there isn't enough room
   */
  ConeVisualization* addCone() {
    if(num_cones < VISUALIZATION_MAX_ITEMS) {
      return &cones[num_cones++];
    }
    return nullptr;
  }

  /*!
   * Add a new path
   * @return A path, or nullptr if there isn't enough room
   */
  PathVisualization* addPath() {
    if(num_paths < VISUALIZATION_MAX_PATHS) {
      auto* path = &paths[num_paths++];
      path->clear();
      return path;
    }
    return nullptr;
  }

  /*!
   * Add a new Mesh
   * @return A mesh, or nullptr if there isn't enough room
   */
   MeshVisualization* addMesh() { 
    if(num_paths < VISUALIZATION_MAX_MESHES) {
      return &meshes[num_meshes++];
    }
    return nullptr;
  }
};

#endif
