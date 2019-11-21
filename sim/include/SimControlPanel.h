/*!
 * @file SimControlPanel.h
 * @brief QT gui for the simulator
 */

#ifndef SIMCONTROLPANEL_H
#define SIMCONTROLPANEL_H

#include <QMainWindow>
#include <thread>
#include "ControlParameters/SimulatorParameters.h"
#include "Graphics3D.h"
#include "RobotInterface.h"
#include "Simulation.h"

#define DEFAULT_TERRAIN_FILE "/default-terrain.yaml"
#define DEFAULT_USER_FILE "/default-user-parameters-file.yaml"

#include <lcm-cpp.hpp>
#include <src/MiniCheetahDebug.h>
#include <leg_control_data_lcmt.hpp>
#include "rs_pointcloud_t.hpp"
#include "heightmap_t.hpp"
#include "traversability_map_t.hpp"
#include "obstacle_visual_t.hpp"
#include "velocity_visual_t.hpp"

namespace Ui {
class SimControlPanel;
}

enum class SimulationWindowState {
  STOPPED,
  RUNNING,
  ERROR
};

class SimControlPanel : public QMainWindow {
  Q_OBJECT

 public:
  explicit SimControlPanel(QWidget* parent = nullptr);
  ~SimControlPanel();



public slots:
  void update_ui();
  void errorCallback(std::string errorMessage);

 private slots:

  void on_startButton_clicked();

  void on_stopButton_clicked();

  void on_joystickButton_clicked();

  void on_driverButton_clicked();

  void on_simulatorTable_cellChanged(int row, int column);

  void on_saveSimulatorButton_clicked();

  void on_loadSimulatorButton_clicked();

  void on_robotTable_cellChanged(int row, int column);

  void on_saveRobotButton_clicked();

  void on_loadRobotButton_clicked();

  void on_goHomeButton_clicked();

  void on_kickButton_clicked();

  void on_userControlTable_cellChanged(int row, int column);

  void on_saveUserButton_clicked();

  void on_loadUserButton_clicked();

  void on_setTerrainButton_clicked();

  void on_hide_floor_checkbox_toggled(bool x) {
    if(_graphicsWindow) {
      _graphicsWindow->setHideFloor(x);
    }
  }

  void on_hide_robot_checkbox_toggled(bool x) {
    if(_graphicsWindow) {
      _graphicsWindow->setHideRobot(x);
    }
  }

  void updateTerrainLabel();

  void loadSimulationParameters(SimulatorControlParameters& params);
  void loadRobotParameters(RobotControlParameters& params);
  void loadUserParameters(ControlParameters& params);

 private:

  std::string getDefaultUserParameterFileName();
  void updateUiEnable();
  bool isStopped() {
    return _state == SimulationWindowState::STOPPED;
  }

  bool isError() {
    return _state == SimulationWindowState::ERROR;
  }

  bool isRunning() {
    return _state == SimulationWindowState::RUNNING;
  }


  std::thread _simThread;
  //bool _started = false;
  SimulationWindowState _state = SimulationWindowState::STOPPED;
  Ui::SimControlPanel* ui;
  Simulation* _simulation = nullptr;
  PeriodicTaskManager* _interfaceTaskManager = nullptr;
  RobotInterface* _robotInterface = nullptr;
  Graphics3D* _graphicsWindow = nullptr;
  SimulatorControlParameters _parameters;
  ControlParameters _userParameters;
  bool _simulationMode = false;
  bool _firstStart = true;
  bool _ignoreTableCallbacks = false;
  bool _loadedUserSettings = false;
  std::string _terrainFileName;

  // Vision data Drawing
  void handleHeightmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, 
      const heightmap_t* msg);
  void heightmapLCMThread() { while (true) { _heightmapLCM.handle(); } }

  void handlePointsLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, 
      const rs_pointcloud_t* msg);
  void pointsLCMThread() { while (true) { _pointsLCM.handle(); } }

  void handleIndexmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, 
      const traversability_map_t* msg);
  void indexmapLCMThread() { while (true) { _indexmapLCM.handle(); } }

  void handleObstacleLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, 
      const obstacle_visual_t* msg);
  void handleVelocityCMDLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, 
      const velocity_visual_t* msg);
  void ctrlVisionLCMThread(){ while(true){ _ctrlVisionLCM.handle();  } }

  void handleSpiDebug(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const leg_control_data_lcmt* msg);

  lcm::LCM _heightmapLCM;
  lcm::LCM _pointsLCM;
  lcm::LCM _indexmapLCM;
  lcm::LCM _ctrlVisionLCM;
  lcm::LCM _miniCheetahDebugLCM;

  std::thread _pointsLCMThread;
  std::thread _heightmapLCMThread;
  std::thread _indexmapLCMThread;
  std::thread _ctrlVisionLCMThread;
  std::thread _miniCheetahDebugLCMThread;

  MiniCheetahDebug _mcDebugWindow;

};

#endif  // SIMCONTROLPANEL_H
