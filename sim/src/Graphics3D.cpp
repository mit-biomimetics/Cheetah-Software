/*! @file Graphics3D.cpp
 *  @brief Visualizer window for simulator
 *
 *  This class displays a window for 3D graphics. It also implements
 * scroll/pan/zoom.
 */

#include "Graphics3D.h"
#include "Utilities/utilities.h"

#include <GL/glut.h>
#include <unistd.h>
#include <iostream>

// background color
static constexpr auto clearColor = windows2000;

// shader program for drawing models where each vertex is assigned a color
static constexpr char vertexShaderSourceColorArray[] = R"(
// inputs:
attribute highp vec3 posAttr;   // position
attribute lowp vec3 colAttr;      // color
attribute highp vec3 normAttr;  // normal
uniform highp mat4 matrix;      // transformation
// outputs:
varying lowp vec4 col;          // color
varying vec3 normal;            // normal
varying vec3 pos_world;         // position
void main() {
  col = vec4(colAttr,1);
  gl_Position = matrix * vec4(posAttr,1);
  normal = (matrix * vec4(normAttr,0)).xyz;
  pos_world = posAttr;
}
)";

// shader program for drawing models where each vertex
static constexpr char vertexShaderSourceSolidColor[] = R"(
// inputs:
attribute highp vec3 posAttr;   // position
uniform lowp vec4 colUniform;          // color
attribute lowp vec3 colAttr;      // color
attribute highp vec3 normAttr;  // normal
uniform highp mat4 matrix;      // transformation
// outputs:
varying lowp vec4 col;          // color
varying vec3 normal;            // normal
varying vec3 pos_world;         // position
void main() {
  col = colUniform;
  gl_Position = matrix * vec4(posAttr,1);
  normal = (matrix * vec4(normAttr,0)).xyz;
  pos_world = posAttr;
}
)";

static constexpr char fragmentShaderSource[] = R"(
varying lowp vec4 col;
varying vec3 pos_world;
varying vec3 normal;
void main() {
  vec3 light_pos = vec3(0,0,4);
  float light_dist = length(light_pos - pos_world);
  vec3 light_color = vec3(1,1,1);
  vec3 mat_ambient = vec3(.1,.1,.1);
  vec3 mat_spec    = vec3(.3,.3,.3);
  vec3 n = normalize(normal);
  vec3 l = normalize(light_pos);
  float angle_thing = clamp( dot(-n,l), 0., 1.);
  gl_FragColor = col + vec4(mat_ambient*light_color*40.*angle_thing/(light_dist*light_dist),0);
}
)";

/*!
 * Initialize a 3D visualization window
 */
Graphics3D::Graphics3D(QWidget *parent)
    : QOpenGLWidget(parent),
      _animating(false),
      _colorArrayProgram(0),
      _frame(0),
      _v0(0, 0, 0),
      _freeCamFilter(1, 60, _v0) {
  std::cout << "[SIM GRAPHICS] New graphics window. \n";

  _r[0] = 0.2422;
  _g[0] = 0.1504;
  _b[0] = 0.6603;
  _r[1] = 0.2810;
  _g[1] = 0.3228;
  _b[1] = 0.9579;
  _r[2] = 0.1786;
  _g[2] = 0.5289;
  _b[2] = 0.9682;
  _r[3] = 0.0689;
  _g[3] = 0.6948;
  _b[3] = 0.8394;
  _r[4] = 0.2161;
  _g[4] = 0.7843;
  _b[4] = 0.5923;
  _r[5] = 0.6720;
  _g[5] = 0.7793;
  _b[5] = 0.2227;
  _r[6] = 0.9970;
  _g[6] = 0.7659;
  _b[6] = 0.2199;
  _r[7] = 0.9769;
  _g[7] = 0.9839;
  _b[7] = 0.0805;
}

Graphics3D::~Graphics3D() {}

/*!
 * Configure the window for displaying cheetah 3
 */
size_t Graphics3D::setupCheetah3(Vec4<float> color, bool useOld) { return _drawList.addCheetah3(color, useOld); }

/*!
 * Configure the window for displaying mini cheetah
 */
size_t Graphics3D::setupMiniCheetah(Vec4<float> color, bool useOld) { return _drawList.addMiniCheetah(color, useOld); }

/*!
 * Update the camera matrix for the current zoom/orbit
 */

void Graphics3D::updateCameraMatrix() {
  _cameraMatrix.setToIdentity();
  _cameraMatrix.perspective(
      60.f, float(size().width()) / float(size().height()), .001f, 50.f);

  if (_arrowsPressed[0]) _ry -= _targetSpeed / 2.f;
  if (_arrowsPressed[1]) _ry += _targetSpeed / 2.f;
  if (_arrowsPressed[2]) _rx += _targetSpeed / 2.f;
  if (_arrowsPressed[3]) _rx -= _targetSpeed / 2.f;
  if (!_rotOrig) {
    _ry = coerce<float>(_ry, -180, 0);
    // velocity in camera coordinates
    // we want the inverse transformation (coordinateRotation goes the opposite
    // way as QMatrix.rotate())
    RotMat<float> R =
        coordinateRotation<float>(CoordinateAxis::Z, deg2rad(_rx)) *
        coordinateRotation<float>(CoordinateAxis::X, deg2rad(_ry));
    Vec3<float> v(_freeCamMove[0], _freeCamMove[1], _freeCamMove[2]);
    v = R * v;

    // integrate and filter
    _freeCamFilter.update(v);
    for (size_t i = 0; i < 3; i++)
      _freeCamPos[i] += _frameTime * _freeCamFilter.get()[i];

    // apply
    _cameraMatrix.rotate(_ry, 1, 0, 0);
    _cameraMatrix.rotate(_rx, 0, 0, 1);
    _cameraMatrix.translate(_freeCamPos[0], _freeCamPos[1], _freeCamPos[2]);
  } else {
    _cameraMatrix.translate(0.f, 0.f, -.45f * _zoom);
    _cameraMatrix.rotate(_ry, 1, 0, 0);
    _cameraMatrix.rotate(_rx, 0, 0, 1);
    _cameraTarget = _cameraTarget * 0.9 + _drawList.getCameraOrigin() * 0.1;
    _cameraMatrix.translate(-_cameraTarget[0],
                            -_cameraTarget[1], 0);
  }
}

void Graphics3D::initializeGL() {
  std::cout << "[Graphics3D] Initialize OpenGL...\n";
  _cameraTarget.setZero();
  initializeOpenGLFunctions();
  // create GPU shaders
  _colorArrayProgram = new QOpenGLShaderProgram(this);
  _colorArrayProgram->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                              vertexShaderSourceColorArray);
  _colorArrayProgram->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                              fragmentShaderSource);
  _colorArrayProgram->link();

  _solidColorProgram = new QOpenGLShaderProgram(this);
  _solidColorProgram->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                              vertexShaderSourceSolidColor);
  _solidColorProgram->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                              fragmentShaderSource);
  _solidColorProgram->link();

  // setup attributes/uniforms (inputs to shaders)
  _posAttrColorArray = (GLuint)_colorArrayProgram->attributeLocation("posAttr");
  _colAttrColorArray = (GLuint)_colorArrayProgram->attributeLocation("colAttr");
  _normAttrColorArray =
      (GLuint)_colorArrayProgram->attributeLocation("normAttr");
  _matrixUniformColorArray =
      (GLuint)_colorArrayProgram->uniformLocation("matrix");

  _posAttrSolidColor = (GLuint)_solidColorProgram->attributeLocation("posAttr");
  _colUniformSolidColor =
      (GLuint)_solidColorProgram->uniformLocation("colUniform");
  _colAttrSolidColor = (GLuint)_solidColorProgram->attributeLocation("colAttr");
  _matrixUniformSolidColor =
      (GLuint)_solidColorProgram->attributeLocation("normAttr");
  _normAttrSolidColor = (GLuint)_solidColorProgram->uniformLocation("matrix");

  // set clear color:
  glClearColor(clearColor[0], clearColor[1], clearColor[2], 0.f);
}

/*-----------------------------------------*
 * Mouse Handlers for Orbit and Zoom       *
 *-----------------------------------------*/

void Graphics3D::mousePressEvent(QMouseEvent *event) {
  _orbiting = true;
  _orbiting_x_start = event->pos().x();
  _orbiting_y_start = event->pos().y();
  _rx_base = _rx;
  _ry_base = _ry;
}

void Graphics3D::mouseMoveEvent(QMouseEvent *event) {
  if (!_orbiting) return;
  _rx = _rx_base + _pixel_to_rad * (event->pos().x() - _orbiting_x_start);
  _ry = _ry_base + _pixel_to_rad * (event->pos().y() - _orbiting_y_start);
}

void Graphics3D::mouseReleaseEvent(QMouseEvent *event) {
  _orbiting = false;
  _rx_base = _rx_base + _pixel_to_rad * (event->pos().x() - _orbiting_x_start);
  _ry_base = _ry_base + _pixel_to_rad * (event->pos().y() - _orbiting_y_start);
}

void Graphics3D::wheelEvent(QWheelEvent *e) {
  if (e->angleDelta().y() > 0) {
    if (_zoom > .1) _zoom = 0.8f * _zoom;
  } else {
    if (_zoom < 100) _zoom = 1.2f * _zoom;
  }
}

void Graphics3D::keyPressEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_Control) {
    _targetSpeed *= .5f;
  } else if (e->key() == Qt::Key_Shift) {
    _targetSpeed *= 2.f;
  }
  if (e->key() == Qt::Key_W)
    _freeCamMove[2] = _targetSpeed;
  else if (e->key() == Qt::Key_S)
    _freeCamMove[2] = -_targetSpeed;

  if (e->key() == Qt::Key_A)
    _freeCamMove[0] = _targetSpeed;
  else if (e->key() == Qt::Key_D)
    _freeCamMove[0] = -_targetSpeed;

  if (e->key() == Qt::Key_R)
    _freeCamMove[1] = -_targetSpeed;
  else if (e->key() == Qt::Key_F)
    _freeCamMove[1] = _targetSpeed;

  if (e->key() == Qt::Key_Up)
    _arrowsPressed[0] = true;
  else if (e->key() == Qt::Key_Down)
    _arrowsPressed[1] = true;

  if (e->key() == Qt::Key_Right)
    _arrowsPressed[2] = true;
  else if (e->key() == Qt::Key_Left)
    _arrowsPressed[3] = true;

  if (e->key() == Qt::Key_V) {
    if (_pause)
      _pause = false;
    else
      _pause = true;
  }

  if(e->key() == Qt::Key_T) {
    _turbo = true;
  }

  if (e->key() == Qt::Key_Tab) {
    _freeCamPos[0] = 0.f;
    _freeCamPos[1] = 0.f;
    _freeCamPos[2] = 0.f;
    _rx = 0.f;
    _ry = 0.f;
  }
}

void Graphics3D::keyReleaseEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_Control) {
    _targetSpeed /= .5f;
  } else if (e->key() == Qt::Key_Shift) {
    _targetSpeed /= 2.f;
  } else if (e->key() == Qt::Key_Space) {
    _rotOrig = !_rotOrig;
  }

  if (e->key() == Qt::Key_W)
    _freeCamMove[2] = 0;
  else if (e->key() == Qt::Key_S)
    _freeCamMove[2] = 0;

  if (e->key() == Qt::Key_A)
    _freeCamMove[0] = 0;
  else if (e->key() == Qt::Key_D)
    _freeCamMove[0] = 0;

  if (e->key() == Qt::Key_R)
    _freeCamMove[1] = 0;
  else if (e->key() == Qt::Key_F)
    _freeCamMove[1] = 0;

  if (e->key() == Qt::Key_Up)
    _arrowsPressed[0] = false;
  else if (e->key() == Qt::Key_Down)
    _arrowsPressed[1] = false;

  if (e->key() == Qt::Key_Right)
    _arrowsPressed[2] = false;
  else if (e->key() == Qt::Key_Left)
    _arrowsPressed[3] = false;

  if(e->key() == Qt::Key_T) {
    _turbo = false;
  }
}

/*!
 * Enable and disable animation
 */
void Graphics3D::setAnimating(bool animating) { _animating = animating; }

void Graphics3D::renderDrawlist() {
  // reload if needed
  if (_drawList.needsReload()) {
    // upload data

    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);

    glGenBuffers(3, _buffID);

    glBindBuffer(GL_ARRAY_BUFFER, _buffID[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _drawList.getSizeOfAllData(),
                 _drawList.getVertexArray(), GL_STATIC_DRAW);
    glVertexAttribPointer(_posAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, _buffID[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _drawList.getSizeOfAllData(),
                 _drawList.getColorArray(), GL_STATIC_DRAW);
    glVertexAttribPointer(_colAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, _buffID[2]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _drawList.getSizeOfAllData(),
                 _drawList.getNormalArray(), GL_STATIC_DRAW);
    glVertexAttribPointer(_normAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    printf("[Graphics 3D] Uploaded data (%f MB)\n",
           _drawList.getGLDataSizeMB());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
  }

  // draw objects with color arrays
  _colorArrayProgram->bind();

  glBindBuffer(GL_ARRAY_BUFFER, _buffID[0]);
  glVertexAttribPointer(_posAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffID[1]);
  glVertexAttribPointer(_colAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffID[2]);
  glVertexAttribPointer(_normAttrColorArray, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  for (size_t id = 0; id < _drawList.getNumObjectsToDraw(); id++) {
    if (_drawList._instanceColor[id].useSolidColor) {
    } else {
      _colorArrayProgram->setUniformValue(
          _matrixUniformColorArray,
          _cameraMatrix * _drawList.getModelKinematicTransform(id) *
              _drawList.getModelBaseTransform(id));

      glDrawArrays(GL_TRIANGLES, _drawList.getGLDrawArrayOffset(id) / 3,
                   _drawList.getGLDrawArraySize(id) / 3);
    }
  }
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  _colorArrayProgram->release();

  // draw objects without color arrays
  _solidColorProgram->bind();

  glBindBuffer(GL_ARRAY_BUFFER, _buffID[0]);
  glVertexAttribPointer(_posAttrSolidColor, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffID[1]);
  glVertexAttribPointer(_colAttrSolidColor, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffID[2]);
  glVertexAttribPointer(_normAttrSolidColor, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);
  for (size_t id = 0; id < _drawList.getNumObjectsToDraw(); id++) {
    if (_drawList._instanceColor[id].useSolidColor) {
      _solidColorProgram->setUniformValue(
          _matrixUniformSolidColor,
          _cameraMatrix * _drawList.getModelKinematicTransform(id) *
              _drawList.getModelBaseTransform(id));
      auto &col = _drawList._instanceColor[id];
      _solidColorProgram->setUniformValue(
          _colUniformSolidColor,
          QVector4D(col.rgba[0], col.rgba[1], col.rgba[2], col.rgba[3]));

      glDrawArrays(GL_TRIANGLES, _drawList.getGLDrawArrayOffset(id) / 3,
                   _drawList.getGLDrawArraySize(id) / 3);
    } else {
    }
  }
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  _solidColorProgram->release();
}

void Graphics3D::configOpenGLPass(int pass) {
  switch (pass) {
    case 0:
      // clear screen
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // update camera math

      glDisable(GL_CULL_FACE);

      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_LESS);
      break;
    case 1:
      glEnable(GL_CULL_FACE);
      glCullFace(GL_FRONT);
      glDepthFunc(GL_ALWAYS);
      break;
    case 2:
      glEnable(GL_CULL_FACE);
      glCullFace(GL_FRONT);
      glDepthFunc(GL_LEQUAL);
      break;
    case 3:
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      glDepthFunc(GL_ALWAYS);
      break;
    case 4:
      glDisable(GL_CULL_FACE);
      glDepthFunc(GL_LEQUAL);
      break;

    default:
      assert(false);
  }
}

void Graphics3D::paintGL() {
  // update joystick:
  _gameController.updateGamepadCommand(_driverCommand);
  if (!_animating) return;
  if (_frame % 60 == 0) {
    qint64 now = QDateTime::currentMSecsSinceEpoch();
    _fps = (60.f * 1000.f / (now - last_frame_ms));
    // std::cout << "FPS: " << _fps << "\n";
    last_frame_ms = now;
  }
  QPainter painter2(this);
  updateCameraMatrix();

  int passes[] = {0, 1, 2, 3, 4};
  for (int pass : passes) {
    // setup pass:
    glShadeModel(GL_SMOOTH);
    configOpenGLPass(pass);

    // do rendering on all passes:
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    renderDrawlist();
    _Additional_Drawing(pass);
    _BoxObstacleDrawing();
  }

  _MeshObstacleDrawing();

  glDisable(GL_DEPTH_TEST);
  painter2.setPen(QColor(100, 100, 100, 255));
  painter2.fillRect(QRect(30, 30, 400, 200), QColor(100, 100, 100, 220));
  QFont font("Monospace", 20);
  painter2.setPen(QColor(210, 100, 100));
  painter2.setFont(font);
  painter2.drawText(QRect(30, 30, 1000, 1000), Qt::AlignLeft,
                    QString(infoString));
  painter2.end();

  ++_frame;
}

void Graphics3D::_MeshObstacleDrawing() {
  glLoadMatrixf(_cameraMatrix.data());
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  glPushMatrix();
  glTranslatef(_drawList.getHeightMapLeftCorner()[0],
               _drawList.getHeightMapLeftCorner()[1],
               _drawList.getHeightMapLeftCorner()[2]);

  double height_min = _drawList.getHeightMapMax();
  double height_max = _drawList.getHeightMapMin();
  double height_gap = height_max - height_min;
  double scaled_height(0.);

  float color_r(0.f);
  float color_g(0.f);
  float color_b(0.f);

  double grid_size(_drawList.getGridSize());
  double height;
  for (int i(0); i < _drawList.getHeightMap().rows(); ++i) {
    glBegin(GL_LINE_STRIP);
    glPushAttrib(GL_COLOR_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (int j(0); j < _drawList.getHeightMap().cols(); ++j) {
      height = _drawList.getHeightMap()(i, j);
      scaled_height = (height - height_min) / height_gap;
      getHeightColor(scaled_height, color_r, color_g, color_b);
      glColor4f(color_r, color_g, color_b, 1.0f);
      glVertex3d(i * grid_size, j * grid_size, height);
    }
    glPopAttrib();
    glDisable(GL_BLEND);
    glEnd();
  }
  for (int j(0); j < _drawList.getHeightMap().cols(); ++j) {
    glBegin(GL_LINE_STRIP);
    glPushAttrib(GL_COLOR_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (int i(0); i < _drawList.getHeightMap().rows(); ++i) {
      height = _drawList.getHeightMap()(i, j);
      scaled_height = (height - height_min) / height_gap;
      getHeightColor(scaled_height, color_r, color_g, color_b);
      glColor4f(color_r, color_g, color_b, 1.0f);
      glVertex3d(i * grid_size, j * grid_size, height);
    }
    glPopAttrib();
    glDisable(GL_BLEND);
    glEnd();
  }
  glPopMatrix();
  glEnable(GL_LIGHTING);
  glPopAttrib();
}

void Graphics3D::getHeightColor(const double &h, float &r, float &g, float &b) {
  // 0 :  0.2422    0.1504    0.6603
  // 1/7:  0.2810    0.3228    0.9579
  // 2/7: 0.1786    0.5289    0.9682
  // 3/7: 0.0689    0.6948    0.8394
  // 4/7: 0.2161    0.7843    0.5923
  // 5/7: 0.6720    0.7793    0.2227
  // 6/7: 0.9970    0.7659    0.2199
  //   1: 0.9769    0.9839    0.0805
  double step(1. / 7.);
  if (h > 6 * step) {
    _SetRGBHeight(h, step, 6, r, g, b);
  } else if (h > 5 * step) {
    _SetRGBHeight(h, step, 5, r, g, b);
  } else if (h > 4 * step) {
    _SetRGBHeight(h, step, 4, r, g, b);
  } else if (h > 3 * step) {
    _SetRGBHeight(h, step, 3, r, g, b);
  } else if (h > 2 * step) {
    _SetRGBHeight(h, step, 2, r, g, b);
  } else if (h > 1 * step) {
    _SetRGBHeight(h, step, 1, r, g, b);
  } else {
    _SetRGBHeight(h, step, 0, r, g, b);
  }
}
void Graphics3D::_SetRGBHeight(const double &h, const double &step,
                               const int &idx, float &r, float &g, float &b) {
  double grade = (h - idx * step) / step;
  r = (_r[idx + 1] - _r[idx]) * grade + _r[idx];
  g = (_g[idx + 1] - _g[idx]) * grade + _g[idx];
  b = (_b[idx + 1] - _b[idx]) * grade + _b[idx];
}

void Graphics3D::_BoxObstacleDrawing() {
  glLoadMatrixf(_cameraMatrix.data());
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  size_t nBox = _drawList.getBoxInfoList().size();
  for (size_t i(0); i < nBox; ++i) {
    glPushMatrix();
    glMultMatrixf(_drawList.getBoxInfoList()[i].frame);  // column major
    _DrawBox(_drawList.getBoxInfoList()[i].depth,
             _drawList.getBoxInfoList()[i].width,
             _drawList.getBoxInfoList()[i].height);
    glPopMatrix();
  }
  glEnable(GL_LIGHTING);
  glPopAttrib();
}

void Graphics3D::_DrawBox(double depth, double width, double height) {
  double x = depth / 2.0;
  double y = width / 2.0;
  double z = height / 2.0;

  glPushAttrib(GL_COLOR_BUFFER_BIT);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glColor4f(_color1[0], _color1[1], _color1[2], 0.7f);

  glBegin(GL_QUADS);
  {
    glVertex3d(x, y, -z);
    glVertex3d(-x, y, -z);
    glVertex3d(-x, y, z);
    glVertex3d(x, y, z);

    glVertex3d(x, -y, -z);
    glVertex3d(x, y, -z);
    glVertex3d(x, y, z);
    glVertex3d(x, -y, z);

    glVertex3d(x, -y, z);
    glVertex3d(-x, -y, z);
    glVertex3d(-x, -y, -z);
    glVertex3d(x, -y, -z);

    glVertex3d(-x, -y, z);
    glVertex3d(-x, y, z);
    glVertex3d(-x, y, -z);
    glVertex3d(-x, -y, -z);

    glVertex3d(-x, -y, -z);
    glVertex3d(-x, y, -z);
    glVertex3d(x, y, -z);
    glVertex3d(x, -y, -z);

    glVertex3d(-x, -y, z);
    glVertex3d(x, -y, z);
    glVertex3d(x, y, z);
    glVertex3d(-x, y, z);
  }  // GL_QUADS
  glEnd();

  glPopAttrib();
  glDisable(GL_BLEND);
}

void Graphics3D::_Additional_Drawing(int pass) {
  glLoadMatrixf(_cameraMatrix.data());
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  _DrawContactForce();
  _DrawContactPoint();

  configOpenGLPass(pass);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  for (size_t i = 0; i < _drawList._visualizationData->num_arrows; i++) {
    _drawArrow(_drawList._visualizationData->arrows[i]);
  }

  for (size_t i = 0; i < _drawList._visualizationData->num_cones; i++) {
    _drawCone(_drawList._visualizationData->cones[i]);
  }

  for (size_t i = 0; i < _drawList._visualizationData->num_blocks; i++) {
    _drawBlock(_drawList._visualizationData->blocks[i]);
  }

  for (size_t i = 0; i < _drawList._visualizationData->num_spheres; i++) {
    _drawSphere(_drawList._visualizationData->spheres[i]);
  }
  glDisable(GL_BLEND);

  for (size_t i = 0; i < _drawList._visualizationData->num_paths; i++) {
    PathVisualization path = _drawList._visualizationData->paths[i];
    glColor4f(path.color[0], path.color[1], path.color[2], path.color[3]);
    glBegin(GL_LINE_STRIP);
    for (size_t j = 0; j < path.num_points; j++) {
      glVertex3d(path.position[j][0], path.position[j][1], path.position[j][2]);
    }
    glEnd();
  }

  glPopAttrib();
  glEnable(GL_LIGHTING);
}

void Graphics3D::_DrawContactForce() {
  glLineWidth(2.0);
  // double scale(0.02);
  double scale(0.0025);

  glPushAttrib(GL_COLOR_BUFFER_BIT);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  for (size_t i(0); i < _drawList.getTotalNumGC(); ++i) {
    glColor4f(0.8f, 0.0f, 0.f, 0.5f);
    Vec3<float> floatForce, floatPos;
    for (size_t j = 0; j < 3; j++) {
      floatPos(j) = _drawList.getGCPos(i)[j];
      floatForce(j) = scale * _drawList.getGCForce(i)[j];
    }
    _drawArrow(floatPos, floatForce, .005, .015, .04);
  }
  glPopAttrib();
  glDisable(GL_BLEND);
}

void Graphics3D::_DrawContactPoint() {
  glPointSize(5);

  glPushAttrib(GL_COLOR_BUFFER_BIT);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  for (size_t i(0); i < _drawList.getTotalNumGC(); ++i) {
    glBegin(GL_POINTS);
    glColor4f(0.8f, 0.0f, 0.1f, 0.3f);

    glVertex3f(_drawList.getGCPos(i)[0], _drawList.getGCPos(i)[1],
               _drawList.getGCPos(i)[2]);

    glEnd();
  }
  glPopAttrib();
  glDisable(GL_BLEND);
}

void Graphics3D::_rotateZtoDirection(const Vec3<float> &direction) {
  float len = direction.norm();
  float dxn = direction(0) / len;
  float dyn = direction(1) / len;
  float dzn = direction(2) / len;

  // Note: We need to create a rotation such that the z-axis points in the
  // direction of the 'direction' argument
  //       Thus, we can realize the rotation with a pure x-y axial rotation
  //       The rotation matrix of the rotated frame 'r' to the current frame 'c'
  //       (c_R_r) has special form. It's 3-rd column in particular has form: [
  //       wy*sin(theta) -wx*sin(theta) (1- cos(theta))]' where theta , [wx wy
  //       0] is the angle, axis of rotation

  // Find the rotation angle by comparing the z-axis of the current and rotated
  // frame
  const double cosTheta = dzn;
  const double theta = acos(cosTheta);
  const double sinTheta = sin(theta);

  // Exploit the special form of the rotation matrix (explained above) for find
  // the axis of rotation
  float rX, rY, rZ;
  if (theta > 0) {
    rX = -dyn / sinTheta;
    rY = dxn / sinTheta;
    rZ = 0;
  } else {
    rX = 0;
    rY = 0;
    rZ = 1;
  }
  glRotatef(ori::rad2deg(theta), rX, rY, rZ);
}

void Graphics3D::_drawSphere(SphereVisualization &sphere) {
  static GLUquadric *quad = gluNewQuadric();
  glPushMatrix();
  _translate(sphere.position);
  _setColor(sphere.color);
  gluSphere(quad, sphere.radius, 16, 16);
  glPopMatrix();
}

void Graphics3D::_drawCone(ConeVisualization &cone) {
  static GLUquadric *quad = gluNewQuadric();
  float len = cone.direction.norm();

  glPushMatrix();
  _setColor(cone.color);
  _translate(cone.point_position);
  _rotateZtoDirection(cone.direction);
  const int detail = 32;
  gluCylinder(quad, 0, cone.radius, len, detail, 1);
  glPopMatrix();
}

void Graphics3D::_drawBlock(BlockVisualization &box) {
  glPushMatrix();
  _translate(box.corner_position);
  _setColor(box.color);
  glRotatef(box.rpy(2), 0, 0, 1);
  glRotatef(box.rpy(1), 0, 1, 0);
  glRotatef(box.rpy(0), 1, 0, 0);

  glScalef(box.dimension(0), box.dimension(1), box.dimension(2));
  glTranslatef(.5, .5, .5);
  _DrawBox(1, 1, 1);
  glPopMatrix();
}

void Graphics3D::_drawArrow(ArrowVisualization &arrow) {
  _setColor(arrow.color);
  _drawArrow(arrow.base_position, arrow.direction, arrow.shaft_width,
             arrow.head_width, arrow.head_length);
}

void Graphics3D::_drawArrow(const Vec3<float> &position,
                            const Vec3<float> &direction, float lineWidth,
                            float headWidth, float headLength) {
  double len = direction.norm();

  glPushMatrix();
  _translate(position);
  _rotateZtoDirection(direction);

  double cylinderLength = len;
  if (cylinderLength > headLength) {
    cylinderLength -= headLength;
  } else {
    headLength = cylinderLength;
    cylinderLength = 0;
  }

  const int detail = 32;
  static GLUquadric *quad = gluNewQuadric();

  // glPolygonMode(GL_FRONT, GL_FILL);
  // Draw Cylinder
  gluCylinder(quad, lineWidth, lineWidth, cylinderLength, detail, 1);

  // Draw Cylinder Base
  glRotated(180, 1, 0, 0);
  gluDisk(quad, 0, lineWidth, detail, detail);
  glRotated(180, 1, 0, 0);

  glTranslatef(0, 0, cylinderLength);
  // Draw Arrowhead
  gluCylinder(quad, headWidth, 0.0f, headLength, detail, detail);

  glRotated(180, 1, 0, 0);
  // Draw Arrowhead Base
  gluDisk(quad, lineWidth, headWidth, detail, detail);
  glPopMatrix();
}
