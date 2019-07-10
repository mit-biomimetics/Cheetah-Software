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
#define DEFAULT_USER_FILE "/default-user.yaml"

namespace Ui {
class SimControlPanel;
}

class SimControlPanel : public QMainWindow {
  Q_OBJECT

 public:
  explicit SimControlPanel(QWidget* parent = nullptr);

  ~SimControlPanel();

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

  void updateTerrainLabel();

  void loadSimulationParameters(SimulatorControlParameters& params);
  void loadRobotParameters(RobotControlParameters& params);
  void loadUserParameters(ControlParameters& params);

 private:
  void updateUiEnable();
  std::thread _simThread;
  bool _started = false;
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
};

#endif  // SIMCONTROLPANEL_H
