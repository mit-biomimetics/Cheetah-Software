/*! @file DrawList.cpp
 *  @brief Data structure to store robot model to be drawn.
 *
 *  Stores all the data (except for joint positions) for the robot.
 *  Knows how to load cheetah robots from file.
 */

#include "DrawList.h"

void DrawList::loadFiles() {
  printf("[DrawList] Load object files...\n");
  std::vector<std::string> names = {
      "c3_body.obj",         "mini_abad.obj",
      "c3_upper_link.obj",   "c3_lower_link.obj",
      "mini_body.obj",       "mini_abad.obj",
      "mini_upper_link.obj", "mini_lower_link.obj",
      "sphere.obj",          "cube.obj"};
  for (const auto& name : names) {
    std::string filename = _baseFileName + name;
    _vertexData.emplace_back();
    _normalData.emplace_back();
    _colorData.emplace_back();
    load_obj_file(filename, _vertexData.back(), _normalData.back());
    if (name == "sphere.obj") {
      setSolidColor(_colorData.back(), _vertexData.back().size(),
                    debugRedColor[0], debugRedColor[1], debugRedColor[2]);
    } else if (name == "cube.obj") {
      setSolidColor(_colorData.back(), _vertexData.back().size(),
                    disgustingGreen[0], disgustingGreen[1], disgustingGreen[2]);
    } else {
      setSolidColor(_colorData.back(), _vertexData.back().size(),
                    defaultRobotColor[0], defaultRobotColor[1],
                    defaultRobotColor[2]);
    }

    _nUnique++;
  }
  _sphereLoadIndex = 8;
  _cubeLoadIndex = 9;
  _miniCheetahLoadIndex = 4;
  _cheetah3LoadIndex = 0;
}
/*!
 * Load the cheetah 3 model and build the draw list.
 * Returns an index number that can later be used to update the position of the
 * robot.
 */
size_t DrawList::addCheetah3(Vec4<float> color, bool useOld, bool canHide) {
  size_t i0 = _cheetah3LoadIndex;
  size_t j0 = _nTotal;

  // set model offsets:
  QMatrix4x4 bodyOffset, abadOffset, lowerOffset, eye;
  eye.setToIdentity();
  QMatrix4x4 upperOffsets[2];

  // body
  bodyOffset.setToIdentity();
  bodyOffset.rotate(90, 1, 0, 0);
  bodyOffset.rotate(90, 0, 0, 1);

  // abad
  abadOffset.setToIdentity();

  // upper link
  upperOffsets[0].setToIdentity();
  upperOffsets[0].rotate(-180, 0, 1, 0);
  upperOffsets[0].translate(0, -.045f, 0);
  upperOffsets[0].rotate(-180, 0, 0, 1);

  upperOffsets[1].setToIdentity();
  upperOffsets[1].rotate(-180, 0, 1, 0);
  upperOffsets[1].translate(0, .045f, 0);

  lowerOffset.setToIdentity();
  lowerOffset.rotate(180, 0, 1, 0);
  lowerOffset.translate(0, 0, 0);

  SolidColor bodyColor, abadColor, link1Color, link2Color;
  bodyColor.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  bodyColor.useSolidColor = true;

  abadColor.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  abadColor.useSolidColor = true;

  link1Color.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  link1Color.useSolidColor = true;

  link2Color.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  link2Color.useSolidColor = true;

  _canBeHidden.push_back(canHide);

  // add bodies
  _objectMap.push_back(i0 + 0);
  _modelOffsets.push_back(bodyOffset);
  _kinematicXform.push_back(eye);
  _instanceColor.push_back(bodyColor);
  _nTotal++;

  for (int i = 0; i < 4; i++) {
    _objectMap.push_back(i0 + 1);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(abadOffset);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(abadColor);

    _objectMap.push_back(i0 + 2);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(upperOffsets[i % 2]);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(link1Color);

    _objectMap.push_back(i0 + 3);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(lowerOffset);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(link2Color);
    _nTotal += 3;
  }

  printf("size of kinematicXform: %lu, j0: %lu\n", _kinematicXform.size(), j0);

  buildDrawList();
  return j0;
}

/*!
 * Load the mini cheetah model and builds the draw list.
 * Returns an index number that can later be used to update the position of the
 * robot.
 * TODO check all this once the mini cheetah dynamics model exists again
 */
size_t DrawList::addMiniCheetah(Vec4<float> color, bool useOld, bool canHide) {
  
  size_t i0 = _miniCheetahLoadIndex;  // todo don't hard code this
  size_t j0 = _nTotal;

  // set model offsets:
  QMatrix4x4 bodyOffset, upper, lower, eye;
  QMatrix4x4 abadOffsets[4];
  eye.setToIdentity();

  // body
  bodyOffset.setToIdentity();

  // abads (todo, check these)
  abadOffsets[0].setToIdentity();  // n
  abadOffsets[0].rotate(-90, 0, 0, 1);
  abadOffsets[0].translate(0, -.0565f, 0);
  abadOffsets[0].rotate(180, 0, 1, 0);

  abadOffsets[1].setToIdentity();  // p
  abadOffsets[1].rotate(-90, 0, 0, 1);
  abadOffsets[1].translate(0, -.0565f, 0);
  abadOffsets[1].rotate(0, 0, 1, 0);

  abadOffsets[2].setToIdentity();  // n
  abadOffsets[2].rotate(90, 0, 0, 1);
  abadOffsets[2].translate(0, -.0565f, 0);
  abadOffsets[2].rotate(0, 0, 1, 0);

  abadOffsets[3].setToIdentity();  // p
  abadOffsets[3].rotate(90, 0, 0, 1);
  abadOffsets[3].translate(0, -.0565f, 0);
  abadOffsets[3].rotate(180, 0, 1, 0);

  // upper
  upper.setToIdentity();
  upper.rotate(-90, 0, 1, 0);

  // lower
  lower.setToIdentity();
  lower.rotate(180, 0, 1, 0);

  SolidColor bodyColor, abadColor, link1Color, link2Color;
  bodyColor.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  bodyColor.useSolidColor = true;

  abadColor.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  abadColor.useSolidColor = true;

  link1Color.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  link1Color.useSolidColor = true;

  link2Color.rgba = useOld ? Vec4<float>(.2, .2, .4, .3) : color;
  link2Color.useSolidColor = true;

  _canBeHidden.push_back(canHide);

  // add objects
  _objectMap.push_back(i0 + 0);
  _modelOffsets.push_back(bodyOffset);
  _kinematicXform.push_back(eye);
  _instanceColor.push_back(bodyColor);
  _nTotal++;

  for (int i = 0; i < 4; i++) {
    _objectMap.push_back(i0 + 1);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(abadOffsets[i]);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(abadColor);

    _objectMap.push_back(i0 + 2);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(upper);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(link1Color);

    _objectMap.push_back(i0 + 3);
    _canBeHidden.push_back(canHide);
    _modelOffsets.push_back(lower);
    _kinematicXform.push_back(eye);
    _instanceColor.push_back(link2Color);
    _nTotal += 3;
  }

  // printf("add mini cheetah (%d) id %ld\n", (int)canHide, j0);
  // for(u32 i = 0; i < _canBeHidden.size(); i++) {
  //   printf(" [%02d] %d\n", i, _canBeHidden[i]);
  // }
  return j0;
}

/*!
 * Adds a checkerboard to the list of drawables.
 * Uses an identity transformation. You must call
 * updateCheckerboardFromCollisionPlane to set the actual transform.
 */
size_t DrawList::addCheckerboard(Checkerboard& checkerBoard, bool scroll) {
  size_t j0 = _nTotal;
  size_t i0 = _nUnique;

  SolidColor checkerColor;
  checkerColor.useSolidColor = false;

  _nUnique++;
  // add the object
  _vertexData.emplace_back();
  _normalData.emplace_back();
  _colorData.emplace_back();
  checkerBoard.computeVertices(_vertexData.back(), _normalData.back(),
                               _colorData.back());
  QMatrix4x4 eye, offset;
  eye.setToIdentity();
  offset.setToIdentity();
  offset.translate(-checkerBoard.getSize()[0] / 2,
                   -checkerBoard.getSize()[1] / 2);
  _modelOffsets.push_back(offset);
  _kinematicXform.push_back(eye);
  _instanceColor.push_back(checkerColor);

  _nTotal++;
  // add the instance
  _objectMap.push_back(i0);
  _canBeHidden.push_back(false);
  if(scroll) {
    _scrollIDs.push_back({j0, checkerBoard.getSize()[0], checkerBoard.getSize()[1]});
  }
  return j0;
}

void DrawList::doScrolling(Vec3<float> cameraPos) {
  for(auto& obj : _scrollIDs) {
    float scrollDiv[2] = {obj.xs/4, obj.ys/4};
    auto& groundXform = getModelKinematicTransform(obj.id);
    groundXform.setToIdentity();
    groundXform.translate( -scrollDiv[0] * (int)(cameraPos[0] / scrollDiv[0]),  -scrollDiv[1] * (int)(cameraPos[1] / scrollDiv[1]));
  }
}

/*!
 * Adds a sphere to the list of drawables.
 */
size_t DrawList::addDebugSphere(float radius) {
  assert(false);
  size_t j0 = _nTotal;

  QMatrix4x4 offset;
  offset.setToIdentity();
  _kinematicXform.push_back(offset);
  offset.scale(radius);
  _modelOffsets.push_back(offset);

  _nTotal++;
  _objectMap.push_back(_sphereLoadIndex);
  _canBeHidden.push_back(false);
  return j0;
}

/*!
 * Rebuilds the drawing list and sets the flag indicating that model data must
 * be reloaded.
 */
void DrawList::buildDrawList() {
  _glVertexData.clear();
  _glColorData.clear();
  _glNormalData.clear();
  _glArrayOffsets.clear();
  _glArraySizes.clear();

  for (size_t i = 0; i < _nUnique; i++) {
    _glArrayOffsets.push_back(_glVertexData.size());
    _glArraySizes.push_back(_vertexData.at(i).size());
    // add the data for the objects
    _glVertexData.insert(_glVertexData.end(), _vertexData.at(i).begin(),
                         _vertexData.at(i).end());
    _glColorData.insert(_glColorData.end(), _colorData.at(i).begin(),
                        _colorData.at(i).end());
    _glNormalData.insert(_glNormalData.end(), _normalData.at(i).begin(),
                         _normalData.at(i).end());
  }

  _reloadNeeded = true;
}

void DrawList::addBox(double depth, double width, double height,
                      const Vec3<double>& pos, const Mat3<double>& ori,
                      bool transparent) {
  if (transparent) {
    BoxInfo tmp;
    tmp.depth = depth;
    tmp.width = width;
    tmp.height = height;

    tmp.frame[3] = 0.;
    tmp.frame[7] = 0.;
    tmp.frame[11] = 0.;
    tmp.frame[15] = 1.;

    for (size_t i(0); i < 3; ++i) {
      for (size_t j(0); j < 3; ++j) {
        tmp.frame[4 * i + j] = ori(j, i);
      }
    }
    for (size_t i(0); i < 3; ++i) tmp.frame[12 + i] = pos[i];

    _box_list.push_back(tmp);
  } else {
    QMatrix4x4 offset;

    // scale box
    offset.setToIdentity();
    offset.scale(depth, width, height);
    _modelOffsets.push_back(offset);

    // move box
    offset.setToIdentity();
    offset.translate(pos[0], pos[1], pos[2]);
    Quat<double> q = ori::rotationMatrixToQuaternion(ori.transpose());
    QQuaternion qq(q[0], q[1], q[2], q[3]);
    offset.rotate(qq);

    _kinematicXform.push_back(offset);

    SolidColor boxColor;
    boxColor.rgba = Vec4<float>(disgustingGreen[0], disgustingGreen[1],
                                disgustingGreen[2], 1.f);
    boxColor.useSolidColor = true;
    _instanceColor.push_back(boxColor);

    _nTotal++;
    _objectMap.push_back(_cubeLoadIndex);
    _canBeHidden.push_back(false);
  }
}

void DrawList::addMesh(double grid_size, const Vec3<double>& left_corner,
                       const DMat<double>& height_map, bool transparent) {
  (void)transparent;

  _grid_size = grid_size;
  _height_map_left_corner = left_corner;
  _height_map = height_map;
  _height_map_min = 1.e5;
  _height_map_max = -1.e5;

  for (int i(0); i < height_map.rows(); ++i) {
    for (int j(0); j < height_map.cols(); ++j) {
      if (height_map(i, j) > _height_map_max) {
        _height_map_max = height_map(i, j);
      }
      if (height_map(i, j) < _height_map_min) {
        _height_map_min = height_map(i, j);
      }
    }
  }
}
