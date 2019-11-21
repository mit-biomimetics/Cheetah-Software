/*! @file obj_loader.cpp
 *  @brief Utility to load .obj files, containing 3D models of robots.
 *
 *  This is now our own version of obj loader
 *  OBJ files should have a triangular mesh, with normals, but no texture
 * coordinates.
 */

#include "obj_loader.h"

#include <assert.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <array>

void load_obj_file(std::string fileName, std::vector<float>& positions,
                   std::vector<float>& normals) {
  std::ifstream file(fileName, std::ios_base::in);
  std::string line;
  std::vector<std::array<float, 3>> verts;
  std::vector<std::array<float, 3>> vns;
  std::vector<std::array<std::array<int, 2>, 3>> faces;

  while (!file.eof()) {
    std::getline(file, line);
    std::istringstream iss(line);
    std::string type;
    iss >> type;
    float v0, v1, v2;

    if (type[0] == '#' || type[0] == 'o') {
      continue;
    } else if (type == "v") {
      iss >> v0;
      iss >> v1;
      iss >> v2;
      verts.push_back({v0, v1, v2});
    } else if (type == "vn") {
      iss >> v0;
      iss >> v1;
      iss >> v2;
      vns.push_back({v0, v1, v2});
    } else if (type == "f") {
      const char* linePtr = line.c_str() + 2;
      int v, n;

      std::array<std::array<int, 2>, 3> face;
      for (int i = 0; i < 3; i++) {
        sscanf(linePtr, "%d//%d", &v, &n);
        face[i][0] = v - 1;
        face[i][1] = n - 1;
        while (linePtr[0] != ' ' && linePtr[0] != 0) ++linePtr;
        while (linePtr[0] == ' ') ++linePtr;
      }
      faces.push_back(face);
    } else if (type == "" || type == "s" || type == "g" || type == "mtllib" || type == "usemtl") {
      // unhandled, but not an error.
    } else {
      printf("ERROR bad obj line %s\n", line.c_str());
    }
  }

  positions.resize(9 * faces.size());
  normals.resize(9 * faces.size());

  uint64_t c = 0;
  for (uint64_t i = 0; i < faces.size(); i++) {
    for (uint64_t j = 0; j < 3; j++) {
      for (uint64_t k = 0; k < 3; k++) {
        positions[c] = verts[faces[i][j][0]][k];
        normals[c++] = vns[faces[i][j][1]][k];
      }
    }
  }
}
