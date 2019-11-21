/*! @file obj_loader.h
 *  @brief Utility to load .obj files, containing 3D models of robots.
 */

#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <string>
#include <vector>

void load_obj_file(std::string fileName, std::vector<float>& positions,
                   std::vector<float>& normals);

#endif  // OBJLOADER_H
