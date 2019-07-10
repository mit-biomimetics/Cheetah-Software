/*! @file DrawList.h
 *  @brief Data structure to store robot model to be drawn.
 *
 *  Stores all the data (except for joint positions) for the robot.
 *  Knows how to load cheetah robots from file.
 *  Will need to support adding a variable number of items (for loading terrain
 * which needs to be drawn)
 *
 *  Also supports having duplicated objects
 */

#ifndef PROJECT_DRAWLIST_H
#define PROJECT_DRAWLIST_H

#include "Checkerboard.h"
#include "Collision/CollisionPlane.h"
#include "Colors.h"
#include "Dynamics/DynamicsSimulator.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/spatial.h"
#include "SimUtilities/VisualizationData.h"
#include "cppTypes.h"
#include "obj_loader.h"
#include "sim_utilities.h"

#include <QMatrix4x4>

#include <stdlib.h>
#include <vector>

class BoxInfo {
 public:
  double depth, width, height;
  float frame[16];  // SE3
                    /*
                       T[0] T[4] T[8]  T[12]
                       T[1] T[5] T[9]  T[13]
                       T[2] T[6] T[10] T[14]
                       T[3] T[7] T[11] T[15] = (0, 0, 0, 1)
                       */
};

struct SolidColor {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vec4<float> rgba;
  bool useSolidColor;
};

class DrawList {
 public:
  VisualizationData *_visualizationData;

  DrawList() {
    _cameraOrigin = Vec3<double>::Zero();
    loadFiles();
  }
  size_t addCheetah3(Vec4<float> color, bool useOld);
  size_t addMiniCheetah(Vec4<float> color, bool useOld);
  void buildDrawList();
  void loadFiles();
  size_t addCheckerboard(Checkerboard &checkerBoard);
  size_t addDebugSphere(float radius);
  void addBox(double depth, double width, double height,
              const Vec3<double> &pos, const Mat3<double> &ori,
              bool transparent);
  void addMesh(double grid_size, const Vec3<double> &left_corner,
               const DMat<double> &height_map, bool transparent);

  /*!
   * Resize to hold size objects.
   */
  void resize(size_t nUniqueObject, size_t nTotalObjects) {
    _nUnique = nUniqueObject;
    _nTotal = nTotalObjects;
    _vertexData.resize(nUniqueObject);
    _normalData.resize(nUniqueObject);
    _colorData.resize(nUniqueObject);
    _offsetXforms.resize(nUniqueObject);
    _objectMap.resize(nTotalObjects);
  }

  /*!
   * Get the total number of objects to be drawn
   */
  size_t getNumObjectsToDraw() { return _nTotal; }

  /*!
   * For the i-th object, get the offset into the model data.
   * For use with the glDrawArrays function.
   * Note that several objects may have the same geometry, so this will return
   * the same value for these objects!
   */
  size_t getGLDrawArrayOffset(size_t i) {
    return _glArrayOffsets.at(_objectMap.at(i));
  }

  /*!
   * For the i-th object, get the size of the model data array
   */
  size_t getGLDrawArraySize(size_t i) {
    return _glArraySizes.at(_objectMap.at(i));
  }

  /*!
   * Get the array containing all vertex data.  Use getGLDrawArrayOffset/Size to
   * get indices and sizes for each object
   */
  float *getVertexArray() { return _glVertexData.data(); }

  /*!
   * Get the array containing all normal data.  Use getGLDrawArrayOffset/Size to
   * get indices and sizes for each object
   */
  float *getNormalArray() { return _glNormalData.data(); }

  size_t getSizeOfAllData() { return _glVertexData.size(); }

  /*!
   * Get the array containing all color data.
   */
  float *getColorArray() { return _glColorData.data(); }

  /*!
   * Get the Qt transformation matrix which should be applied to the model
   * geometry This is to correct for errors when exporting parts and shifts them
   * all to the origin.
   */
  QMatrix4x4 &getModelBaseTransform(size_t i) { return _modelOffsets[i]; }

  /*!
   * Get the Qt transformation matrix which should be applied to move a model
   * from the origin to where it should be in the world
   */
  QMatrix4x4 &getModelKinematicTransform(size_t i) {
    return _kinematicXform[i];
  }

  /*!
   * Get size of data used by the GPU in megabytes
   * For debugging
   */
  float getGLDataSizeMB() {
    size_t bytes =
        _glColorData.size() + _glNormalData.size() + _glVertexData.size();
    bytes = bytes * sizeof(float);
    return (float)bytes / float(1 << 20);
  }

  /*!
   * Returns true a single time if we have changed geometries and need to
   * reload.
   */
  bool needsReload() {
    if (_reloadNeeded) {
      _reloadNeeded = false;
      return true;
    }
    return false;
  }

  /*!
   * Update the position of a robot's bodies using the result of a dynamics
   * simulation. Doesn't run the simulator - just pulls _Xa from the
   * DynamicsSimulator
   * @param model  : the simulator
   * @param id     : the id returned from the loadCheetah3 or loadMiniCheetah
   * function.
   */
  template <typename T>
  void updateRobotFromModel(DynamicsSimulator<T> &model, size_t id,
                            bool updateOrigin = false) {
    for (size_t modelID = 5, graphicsID = id; modelID < model.getNumBodies();
         modelID++, graphicsID++) {
      _kinematicXform.at(graphicsID) =
          spatialTransformToQT(model.getModel()._Xa.at(modelID));
    }

    if (updateOrigin) {
      _cameraOrigin = model.getState().bodyPosition.template cast<T>();
    }
  }

  /*!
   * Update the additional information drawn by GUI
   * Doesn't run the simulator
   *  - just pulls contact (or other in the future) data from the
   * DynamicsSimulator
   * @param model  : the simulator
   */
  template <typename T>
  void updateAdditionalInfo(DynamicsSimulator<T> &model) {
    if (_additionalInfoFirstVisit) {
      _nTotalGC = model.getTotalNumGC();
      _cp_touch.resize(_nTotalGC, false);
      _cp_pos.resize(_nTotalGC);
      _cp_force.resize(_nTotalGC);
      std::vector<double> tmp(3);
      for (size_t i(0); i < _nTotalGC; ++i) {
        _cp_pos[i] = tmp;
        _cp_force[i] = tmp;
      }
      _additionalInfoFirstVisit = false;
    }

    for (size_t i(0); i < _nTotalGC; ++i) {
      // TODO: check touch boolean
      _cp_touch[i] = true;
      for (size_t j(0); j < 3; ++j) {
        _cp_pos[i][j] = model.getModel()._pGC[i][j];
        _cp_force[i][j] = model.getContactForce(i)[j];
      }
    }
  }

  /*!
   * Updates the position of a checkerboard to match an infinite collision plane
   * The infinite collision plane only specifies orientation, so we
   * @param model : the collision plane
   * @param id    : the id retured from creating the checkerboard
   */
  template <typename T>
  void updateCheckerboardFromCollisionPlane(CollisionPlane<T> &model,
                                            size_t id) {
    // Mat4<T> H = sxformToHomogeneous(model.getLocation());
    _kinematicXform.at(id) = spatialTransformToQT(model.getLocation());
  }

  template <typename T>
  void updateCheckerboard(T height, size_t id) {
    QMatrix4x4 H;
    H.setToIdentity();
    H.translate(0., 0., height);
    _kinematicXform.at(id) = H;
  }

  template <typename T>
  void updateDebugSphereLocation(Vec3<T> &position, size_t id) {
    QMatrix4x4 H;
    H.setToIdentity();
    H.translate(position[0], position[1], position[2]);
    _kinematicXform.at(id) = H;
  }

  /*!
   * Fill color data with a solid color
   */
  static void setSolidColor(std::vector<float> &data, size_t size, float r,
                            float g, float b) {
    data.clear();
    data.resize(size);

    if ((size % 3) != 0) {
      throw std::runtime_error("setSolidColor invalid size");
    }

    for (size_t i = 0; i < size / 3; i++) {
      data[i * 3 + 0] = r;
      data[i * 3 + 1] = g;
      data[i * 3 + 2] = b;
    }
  }

  /* Get Functions */
  const size_t &getTotalNumGC() { return _nTotalGC; }
  const std::vector<double> &getGCPos(size_t idx) { return _cp_pos[idx]; }
  const std::vector<double> &getGCForce(size_t idx) { return _cp_force[idx]; }
  const std::vector<BoxInfo> &getBoxInfoList() { return _box_list; }

  const DMat<double> &getHeightMap() { return _height_map; }
  const Vec3<double> &getHeightMapLeftCorner() {
    return _height_map_left_corner;
  }
  const double &getHeightMapMax() { return _height_map_max; }
  const double &getHeightMapMin() { return _height_map_min; }
  const double &getGridSize() { return _grid_size; }

  const Vec3<double> &getCameraOrigin() { return _cameraOrigin; }
  vectorAligned<SolidColor> _instanceColor;
  std::vector<QMatrix4x4> _kinematicXform;

 private:
  size_t _nUnique = 0, _nTotal = 0;
  std::vector<std::vector<float>> _vertexData;
  std::vector<std::vector<float>> _normalData;
  std::vector<std::vector<float>> _colorData;
  vectorAligned<Mat4<float>>
      _offsetXforms;  // these are NOT coordinate transformations!
  std::string _baseFileName = "../resources/";

  std::vector<size_t> _objectMap;

  std::vector<size_t> _glArrayOffsets;
  std::vector<size_t> _glArraySizes;

  std::vector<float> _glVertexData;
  std::vector<float> _glNormalData;
  std::vector<float> _glColorData;

  std::vector<QMatrix4x4> _modelOffsets;

  bool _reloadNeeded = false;
  bool _additionalInfoFirstVisit = true;

  size_t _nTotalGC = 0;
  std::vector<bool> _cp_touch;
  std::vector<std::vector<double>> _cp_pos;
  std::vector<std::vector<double>> _cp_force;
  std::vector<BoxInfo> _box_list;

  double _grid_size;
  Vec3<double> _height_map_left_corner;
  DMat<double> _height_map;
  double _height_map_max, _height_map_min;

  Vec3<double> _cameraOrigin;

  size_t _cheetah3LoadIndex = 0, _miniCheetahLoadIndex = 0,
         _sphereLoadIndex = 0, _cubeLoadIndex = 0;
};

#endif  // PROJECT_DRAWLIST_H
