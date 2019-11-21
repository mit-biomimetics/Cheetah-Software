/*! @file Checkerboard.cpp
 *  @brief 3D Plane with Checkerboard Pattern for graphics
 *
 */

#include "Checkerboard.h"

#include <stdio.h>

/*!
 * Construct a new checkerboard at the origin
 * @param xSize : width, meters
 * @param ySize : height, meters
 * @param xSquares : number of squares along x (must be at least 1)
 * @param ySquares : number of squares along y (must be at least 2)
 */
Checkerboard::Checkerboard(float xSize, float ySize, size_t xSquares,
                           size_t ySquares) {
  _size[0] = xSize;
  _size[1] = ySize;
  _squares[0] = xSquares;
  _squares[1] = ySquares;
  //  computeVertices();
}

/*!
 * Update the checkerboard vertices for the current sizees
 */
void Checkerboard::computeVertices(std::vector<float>& vertices,
                                   std::vector<float>& normals,
                                   std::vector<float>& colors) {
  // 3 channels, n*m squares, 2 triangles, 3 vertices per triangle
  size_t dataSize = 3 * _squares[0] * _squares[1] * 3 * 2;
  vertices.clear();
  normals.clear();
  colors.clear();
  vertices.reserve(dataSize);
  normals.reserve(dataSize);
  colors.reserve(dataSize);

  const float xStep = _size[0] / (float)_squares[0];
  const float yStep = _size[1] / (float)_squares[1];

  for (size_t i = 0; i < _squares[0]; i++) {
    for (size_t j = 0; j < _squares[1]; j++) {
      // add normals and color
      const float* squareColor = (i + j) % 2 ? _darkColor : _lightColor;
      for (size_t k = 0; k < 6; k++) {
        for (size_t m = 0; m < 3; m++) {
          colors.push_back(squareColor[m]);
          normals.push_back(m == 2 ? 1 : 0);
        }
      }

      const float xStart = xStep * i;
      const float yStart = yStep * j;

      // tri 1
      vertices.push_back(xStart);
      vertices.push_back(yStart);
      vertices.push_back(0);

      vertices.push_back(xStart);
      vertices.push_back(yStart + yStep);
      vertices.push_back(0);

      vertices.push_back(xStart + xStep);
      vertices.push_back(yStart);
      vertices.push_back(0);

      // tri 2
      vertices.push_back(xStart);
      vertices.push_back(yStart + yStep);
      vertices.push_back(0);

      vertices.push_back(xStart + xStep);
      vertices.push_back(yStart + yStep);
      vertices.push_back(0);

      vertices.push_back(xStart + xStep);
      vertices.push_back(yStart);
      vertices.push_back(0);
    }
  }
}