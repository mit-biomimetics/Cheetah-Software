/*! @file Graphics3D.h
 *  @brief Visualizer window for simulator
 *
 *  This class displays a window for 3D graphics. It also implements
 * scroll/pan/zoom. This uses OpenGL and QT.
 */

#ifndef PROJECT_GRAPHICS3D_H
#define PROJECT_GRAPHICS3D_H

#include "DrawList.h"
#include "Math/FirstOrderIIRFilter.h"
#include "obj_loader.h"

#include <QMatrix4x4>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLPaintDevice>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <QPainter>
#include <QWindow>
#include "SimUtilities/VisualizationData.h"

#include <QDateTime>
#include <QMouseEvent>
#include <QScreen>
#include <QWheelEvent>

#include <mutex>
#include "GameController.h"

class Graphics3D : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT
 friend class SimControlPanel;
 public:
  explicit Graphics3D(QWidget *parent = 0);
  virtual ~Graphics3D();

  void setAnimating(bool animating);
  size_t setupCheetah3(Vec4<float> color, bool useOld, bool canHide);
  size_t setupMiniCheetah(Vec4<float> color, bool useOld, bool canHide);

  void lockGfxMutex() { _gfxMutex.lock(); }

  void unlockGfxMutex() { _gfxMutex.unlock(); }

  // set robot state
  double _fps = 0;
  DrawList _drawList;
  char infoString[200] = "";

  GamepadCommand &getDriverCommand() { return _driverCommand; }
  GameController &getGameController() { return _gameController; }

  void resetGameController() { _gameController.findNewController(); }

  bool IsPaused() { return _pause; }
  bool wantTurbo() { return _turbo; }
  bool wantSloMo() { return _sloMo; }


  void setHideFloor(bool x);
  void setHideRobot(bool x);

 protected:
  void initializeGL() override;
  // void resizeGL(int w, int h) override;
  void paintGL() override;
  // bool event(QEvent *event) override;

  // void exposeEvent(QExposeEvent *event) override; ??

  // mouse callbacks for orbit and zoom
  void mousePressEvent(QMouseEvent *event) override;
  void mouseMoveEvent(QMouseEvent *event) override;
  void mouseReleaseEvent(QMouseEvent *event) override;
  void wheelEvent(QWheelEvent *e) override;
  void keyReleaseEvent(QKeyEvent *e) override;
  void keyPressEvent(QKeyEvent *event) override;


  float _color1[3] = {0.364784, 0.513401, 0.952230};
  float _color2[3] = {0.553970, 0.477397, 0.628871};
  float _color3[3] = {0.335223, 0.768230, 0.277775};

  bool show_floor = true;
  bool show_robot = true;

 private:
  GameController _gameController;
  GamepadCommand _driverCommand;

  std::mutex _gfxMutex;
  void scrollGround();
  void updateCameraMatrix();
  void renderDrawlist();
  void configOpenGLPass(int pass);
  void _BoxObstacleDrawing();
  void _MeshObstacleDrawing();
  void _DrawBox(double depth, double width, double height);
  void _Additional_Drawing(int pass);
  void _DrawContactForce();
  void _DrawContactPoint();

  void _drawArrow(ArrowVisualization &arrow);
  void _drawBlock(BlockVisualization &box);
  void _drawSphere(SphereVisualization &sphere);
  void _drawCone(ConeVisualization &cone);
  void _drawMesh(MeshVisualization &mesh);
  
  void _rotateZtoDirection(const Vec3<float> &direction);
  void _setColor(const Vec4<float> &color) {
    glColor4f(color(0), color(1), color(2), color(3));
  }
  void _translate(const Vec3<float> &position) {
    glTranslatef(position(0), position(1), position(2));
  }
  void _drawArrow(const Vec3<float> &base, const Vec3<float> &direction,
                  float lineWidth, float headWidth, float headLength);
  bool _animating;

  // attributes for shader program
  GLuint _posAttrColorArray;        // position of vertex
  GLuint _colAttrColorArray;        // color of vertex
  GLuint _matrixUniformColorArray;  // transformation matrix
  GLuint _normAttrColorArray;       // vertex normal

  GLuint _posAttrSolidColor;     // position of vertex
  GLuint _colUniformSolidColor;  // color of vertex
  GLuint _colAttrSolidColor;
  GLuint _matrixUniformSolidColor;  // transformation matrix
  GLuint _normAttrSolidColor;       // vertex normal

  GLuint _buffID[3];

  // shader programs
  QOpenGLShaderProgram *_colorArrayProgram;
  QOpenGLShaderProgram *_solidColorProgram;

  // frame count
  int _frame;
  // time of last frame
  qint64 last_frame_ms = 0;

  // UI orbit/zoom variables
  bool _orbiting = false;
  int _orbiting_x_start = 0;
  int _orbiting_y_start = 0;
  float _rx_base = 0;
  float _ry_base = 0;
  float _rx = 0;
  float _ry = -34;
  float _pixel_to_rad = .3f;
  float _zoom = 3.0;

  bool _rotOrig = true;
  bool _turbo = false;
  bool _sloMo = false;

  QMatrix4x4 _cameraMatrix;
  Vec3<float> _v0;
  Vec3<double> _cameraTarget;
  FirstOrderIIRFilter<Vec3<float>, float> _freeCamFilter;

  float _freeCamMove[3] = {0, 0, 0};
  float _freeCamPos[3] = {0.f, 0.f, 0.f};
  float _frameTime = 1.f / 60.f;

  bool _arrowsPressed[4] = {false, false, false, false};

  float _targetSpeed = 2;

  float _r[8];
  float _g[8];
  float _b[8];

  void getHeightColor(const double &height, float &r, float &g, float &b);
  void _SetRGBHeight(const double &h, const double &step, const int &idx,
                     float &r, float &g, float &b);

  bool _pause = false;

  // Vision data visualization
  Vec3<float> _points[5001];
  size_t _num_points = 5001;
  bool _pointcloud_data_update = false;

  size_t x_size = 100;
  size_t y_size = 100;
  DMat<float> _map;
  DMat<int> _idx_map;
  Vec3<float> _pos;

  bool _heightmap_data_update = false;
  bool _indexmap_data_update = false;
  
  Vec3<float> _vel_cmd_dir, _vel_cmd_pos;
  bool _vel_cmd_update = false;

  vectorAligned< Vec3<double> > _obs_list;
  float _obs_sigma = 0.15;
  float _obs_height;
  bool _obstacle_update = false;

  void _drawHeightMap();
  void _drawVelArrow();
  void _drawObstacleField();
};

#endif  // PROJECT_GRAPHICS3D_H
