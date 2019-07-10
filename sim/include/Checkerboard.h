/*! @file Checkerboard.cpp
 *  @brief 3D Plane with Checkerboard Pattern
 *
 * Oriented with the center at the origin and a 0,0,1 normal
 */

#ifndef PROJECT_CHECKERBOARD_H
#define PROJECT_CHECKERBOARD_H

#include "Colors.h"
#include "cTypes.h"

#include <vector>

class Checkerboard {
 public:
  Checkerboard(float xSize, float ySize, size_t xSquares, size_t ySquares);

  /*!
   * Set the color of the dark squares
   * Do this before computeVertices
   */
  void setDarkColor(const float* dark) { _darkColor = dark; }

  /*!
   * Set the color of the light squares
   * Do this before computeVertices
   */
  void setLightColor(const float* light) { _lightColor = light; }

  /*!
   * Get an array of x,y size (in the plane's coordinates, meters)
   */
  const float* getSize() { return _size; }

  void computeVertices(std::vector<float>& vertices,
                       std::vector<float>& normals, std::vector<float>& colors);

 private:
  const float* _darkColor = checkerboardDark;
  const float* _lightColor = checkerboardLight;
  float _size[2];
  size_t _squares[2];
};

#endif  // PROJECT_CHECKERBOARD_H
