#include "SimControlPanel.h"
#include <ControlParameters/ControlParameters.h>
#include <QFileDialog>
#include <QMessageBox>
#include "ui_SimControlPanel.h"


/*!
 * Display an error messagebox with the given text
 */
static void createErrorMessage(const std::string& text) {
  QMessageBox mb;
  mb.setText(QString(text.c_str()));
  mb.exec();
}

/*!
 * Display control parameters in a qtable.
 */
static void updateQtableWithParameters(ControlParameters& params,
                                       QTableWidget& table) {
  table.setRowCount((s32)params.collection._map.size());
  table.setColumnCount(2);

  s32 i = 0;
  for (auto& kv : params.collection._map) {
    (void)kv;
    for (s32 col = 0; col < 2; col++) {
      QTableWidgetItem* cell = table.item(i, col);
      if (!cell) {
        cell = new QTableWidgetItem;
        table.setItem(i, col, cell);
      }
    }

    table.item(i, 0)->setText(QString(kv.first.c_str()));
    table.item(i, 1)->setText(QString(kv.second->toString().c_str()));
    i++;
  }
}

/*!
 * Init sim window
 */
SimControlPanel::SimControlPanel(QWidget* parent)
    : QMainWindow(parent),
      ui(new Ui::SimControlPanel),
      _userParameters("user-parameters"),
      _terrainFileName(getConfigDirectoryPath() + DEFAULT_TERRAIN_FILE) {

  ui->setupUi(this); // QT setup
  updateUiEnable();  // enable/disable buttons as needed.
  updateTerrainLabel(); // display name of loaded terrain file

  // attempt to load default user settings.
  _loadedUserSettings = true;

  try {
    _userParameters.defineAndInitializeFromYamlFile(getConfigDirectoryPath() + DEFAULT_USER_FILE);
  } catch (std::runtime_error& ex) {
    _loadedUserSettings = false;
  }

  if(!_loadedUserSettings) {
    printf("[SimControlPanel] Failed to load default user settings!\n");
  } else {
    // display user settings in qtable if we loaded successfully
    loadUserParameters(_userParameters);
  }

  // load simulator parameters
  printf("[SimControlPanel] Init simulator parameters...\n");
  _parameters.initializeFromYamlFile(getConfigDirectoryPath() +
                                     SIMULATOR_DEFAULT_PARAMETERS);
  if (!_parameters.isFullyInitialized()) {
    printf(
        "[ERROR] Simulator parameters are not fully initialized.  You forgot: "
        "\n%s\n",
        _parameters.generateUnitializedList().c_str());
    throw std::runtime_error("simulator not initialized");
  } else {
    printf("\tsim parameters are all good\n");
  }
  loadSimulationParameters(_parameters);
}

SimControlPanel::~SimControlPanel() {
  delete _simulation;
  delete _interfaceTaskManager;
  delete _robotInterface;
  delete _graphicsWindow;
  delete ui;
}

/*!
 * Enable/disable buttons as needed based on what is running
 */
void SimControlPanel::updateUiEnable() {
  ui->startButton->setEnabled(!_started);
  ui->stopButton->setEnabled(_started);
  ui->joystickButton->setEnabled(_started);
  ui->robotTable->setEnabled(_started);
  ui->goHomeButton->setEnabled(_started);
}

/*!
 * Update the name of the loaded terrain file label
 */
void SimControlPanel::updateTerrainLabel() {
  ui->terrainFileLabel->setText(QString(_terrainFileName.c_str()));
}

/*!
 * Start a simulation/robot run
 */
void SimControlPanel::on_startButton_clicked() {
  // get robot type
  RobotType robotType;

  if (ui->cheetah3Button->isChecked()) {
    robotType = RobotType::CHEETAH_3;
  } else if (ui->miniCheetahButton->isChecked()) {
    robotType = RobotType::MINI_CHEETAH;
  } else {
    createErrorMessage("Error: you must select a robot");
    return;
  }

  // get run type
  if (!ui->simulatorButton->isChecked() && !ui->robotButton->isChecked()) {
    createErrorMessage(
        "Error: you must select either robot or simulation mode");
    return;
  }

  _simulationMode = ui->simulatorButton->isChecked();

  // graphics
  printf("[SimControlPanel] Initialize Graphics...\n");
  _graphicsWindow = new Graphics3D();
  _graphicsWindow->show();
  _graphicsWindow->resize(1280, 720);

  if (_simulationMode) {
    // run a simulation
    printf("[SimControlPanel] Initialize simulator...\n");
    _simulation = new Simulation(robotType, _graphicsWindow, _parameters, _userParameters);
    loadSimulationParameters(_simulation->getSimParams());
    loadRobotParameters(_simulation->getRobotParams());

    // terrain
    printf("[SimControlParameter] Load terrain...\n");
    _simulation->loadTerrainFile(_terrainFileName);

    // start sim
    _simThread = std::thread([this]() { _simulation->runAtSpeed(); });

    // graphics start
    _graphicsWindow->setAnimating(true);
  } else {
    printf("[SimControlPanel] Init Robot Interface...\n");
    _interfaceTaskManager = new PeriodicTaskManager;
    _robotInterface =
        new RobotInterface(robotType, _graphicsWindow, _interfaceTaskManager, _userParameters);
    loadRobotParameters(_robotInterface->getParams());
    _robotInterface->startInterface();
    _graphicsWindow->setAnimating(true);
  }

  _started = true;
  updateUiEnable();
}

/*!
 * Stop the currently running simulation or robot connection
 */
void SimControlPanel::on_stopButton_clicked() {
  if (_simulation) {
    _simulation->stop();
    _simThread.join();
  } else {
    _robotInterface->stopInterface();
  }

  if (_graphicsWindow) {
    _graphicsWindow->setAnimating(false);
    _graphicsWindow->hide();
  }

  delete _interfaceTaskManager;
  delete _robotInterface;
  delete _simulation;
  delete _graphicsWindow;

  _simulation = nullptr;
  _graphicsWindow = nullptr;
  _robotInterface = nullptr;
  _interfaceTaskManager = nullptr;

  _started = false;
  updateUiEnable();
}



/*!
 * Populate the simulator qtable parameters
 */
void SimControlPanel::loadSimulationParameters(
    SimulatorControlParameters& params) {
  _ignoreTableCallbacks = true;
  updateQtableWithParameters(params, *ui->simulatorTable);
  _ignoreTableCallbacks = false;
}

/*!
 * Populate the robot qtable parameters
 */
void SimControlPanel::loadRobotParameters(RobotControlParameters& params) {
  _ignoreTableCallbacks = true;
  updateQtableWithParameters(params, *ui->robotTable);
  _ignoreTableCallbacks = false;
}

/*!
 * Populate the robot qtable parameters
 */
void SimControlPanel::loadUserParameters(ControlParameters& params) {
  _ignoreTableCallbacks = true;
  updateQtableWithParameters(params, *ui->userControlTable);
  _ignoreTableCallbacks = false;
}

/*!
 * Attempt to reset the joystick if a new one is connected
 */
void SimControlPanel::on_joystickButton_clicked() {
  _graphicsWindow->resetGameController();
}

void SimControlPanel::on_driverButton_clicked() {}

/*!
 * Respond to a change in the simulator table.
 */
void SimControlPanel::on_simulatorTable_cellChanged(int row, int column) {
  if (_ignoreTableCallbacks) return;

  // we only allow values to change, which are in column 1
  if (column != 1) {
    return;
  }

  // get the name of the parameter....
  auto cell = ui->simulatorTable->item(row, 0);
  std::string cellName = cell->text().toStdString();

  if (cellName == "") {
    return;
  }

  // get the parameters
  auto& parameter = _parameters.collection.lookup(cellName);
  ControlParameterValueKind kind = parameter._kind;
  ControlParameterValue oldValue = parameter.get(kind);

  bool success = true;

  // attempt to set, based on string.
  try {
    parameter.setFromString(
        ui->simulatorTable->item(row, 1)->text().toStdString());
  } catch (std::exception& e) {
    success = false;
  }

  // if it fails (bad user input string), restore to the old value
  if (!success) {
    printf("[ERROR] invalid data, restoring old data!\n");
    // set parameter value
    parameter.set(oldValue, kind);

    assert(!_ignoreTableCallbacks);

    // manually fix the table
    _ignoreTableCallbacks = true;
    ui->simulatorTable->item(row, 1)->setText(
        QString(_parameters.collection.lookup(cellName).toString().c_str()));
    _ignoreTableCallbacks = false;
  } else {
    // this update "rewrites" the value in the table.  If it's an integer, it kills any
    // decimal.  If it's a float, it puts in scientific notation if needed.
    _ignoreTableCallbacks = true;
    ui->simulatorTable->item(row, 1)->setText(
        QString(parameter.toString().c_str()));
    _ignoreTableCallbacks = false;
  }
}

/*!
 * Save simulation config to file
 */
void SimControlPanel::on_saveSimulatorButton_clicked() {
  QString fileName = QFileDialog::getSaveFileName(
      nullptr, ("Save Simulator Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  _parameters.lockMutex();
  _parameters.writeToYamlFile(fileName.toStdString());
  _parameters.unlockMutex();
}

/*!
 * Load simulation config from file
 */
void SimControlPanel::on_loadSimulatorButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(
      nullptr, ("Load Simulator Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  };
  _parameters.collection.clearAllSet();
  _parameters.initializeFromYamlFile(fileName.toStdString());
  if (!_parameters.collection.checkIfAllSet()) {
    printf(
        "new settings file %s doesn't contain the following simulator "
        "parameters:\n%s\n",
        fileName.toStdString().c_str(),
        _parameters.generateUnitializedList().c_str());
    throw std::runtime_error("bad new settings file");
  }
  loadSimulationParameters(_parameters);
  _parameters.unlockMutex();
}

void SimControlPanel::on_robotTable_cellChanged(int row, int column) {
  if (_ignoreTableCallbacks) return;
  if (column != 1) {
    return;
  }

  auto cell = ui->robotTable->item(row, 0);
  std::string cellName = cell->text().toStdString();

  if (cellName == "") {
    return;
  }

  auto& parameter = (_simulationMode ? _simulation->getRobotParams()
                                     : _robotInterface->getParams())
                        .collection.lookup(cellName);
  ControlParameterValueKind kind = parameter._kind;
  ControlParameterValue oldValue = parameter.get(kind);

  bool success = true;

  try {
    parameter.setFromString(ui->robotTable->item(row, 1)->text().toStdString());
  } catch (std::exception& e) {
    success = false;
  }

  if (!success) {
    printf("[ERROR] invalid data, restoring old data!\n");
    parameter.set(oldValue, kind);

    assert(!_ignoreTableCallbacks);

    _ignoreTableCallbacks = true;
    ui->robotTable->item(row, 1)->setText(
        QString(parameter
                    .toString()
                    .c_str()));
    _ignoreTableCallbacks = false;
  } else {
    if (_simulationMode) {
      if (_simulation->isRobotConnected()) {
        _simulation->sendControlParameter(
            cellName, parameter.get(parameter._kind), parameter._kind, false);
      }

      _ignoreTableCallbacks = true;
      ui->robotTable->item(row, 1)->setText(
          QString(parameter.toString().c_str()));
      _ignoreTableCallbacks = false;
    } else {
      _robotInterface->sendControlParameter(
          cellName, parameter.get(parameter._kind), parameter._kind, false);
    }
  }
}

void SimControlPanel::on_saveRobotButton_clicked() {
  QString fileName = QFileDialog::getSaveFileName(
      nullptr, ("Save Robot Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }
  _simulation->getRobotParams().writeToYamlFile(fileName.toStdString());
}

void SimControlPanel::on_loadRobotButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(
      nullptr, ("Load Quadruped Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  if (_simulationMode) {
    _simulation->getRobotParams().lockMutex();
    _simulation->getRobotParams().collection.clearAllSet();
    _simulation->getRobotParams().initializeFromYamlFile(
        fileName.toStdString());
    if (!_simulation->getRobotParams().collection.checkIfAllSet()) {
      printf(
          "new settings file %s doesn't contain the following robot "
          "parameters:\n%s\n",
          fileName.toStdString().c_str(),
          _simulation->getRobotParams().generateUnitializedList().c_str());
      throw std::runtime_error("bad new settings file");
    }
    loadRobotParameters(_simulation->getRobotParams());

    if (_simulation->isRobotConnected()) {
      for (auto& kv : _simulation->getRobotParams().collection._map) {
        _simulation->sendControlParameter(
            kv.first, kv.second->get(kv.second->_kind), kv.second->_kind, false);
      }
    }
    _simulation->getRobotParams().unlockMutex();
  } else {
    _robotInterface->getParams().lockMutex();
    _robotInterface->getParams().collection.clearAllSet();
    _robotInterface->getParams().initializeFromYamlFile(fileName.toStdString());
    if (!_robotInterface->getParams().collection.checkIfAllSet()) {
      printf(
          "new settings file %s doesn't contain the following robot "
          "parameters:\n%s\n",
          fileName.toStdString().c_str(),
          _robotInterface->getParams().generateUnitializedList().c_str());
      throw std::runtime_error("bad new settings file");
    }
    loadRobotParameters(_robotInterface->getParams());

    for (auto& kv : _robotInterface->getParams().collection._map) {
      _robotInterface->sendControlParameter(
          kv.first, kv.second->get(kv.second->_kind), kv.second->_kind, false);
    }

    _robotInterface->getParams().unlockMutex();
  }
}

void SimControlPanel::on_setTerrainButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(
      nullptr, ("Load Terrain Definition"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  _terrainFileName = fileName.toStdString();
  updateTerrainLabel();
}

void SimControlPanel::on_userControlTable_cellChanged(int row, int column) {
  if (_ignoreTableCallbacks) return;
  if (column != 1) {
    return;
  }

  auto cell = ui->userControlTable->item(row, 0);
  std::string cellName = cell->text().toStdString();

  if (cellName == "") {
    return;
  }

  auto& parameter = _userParameters.collection.lookup(cellName);
//  auto& parameter = (_simulationMode ? _simulation->getRobotParams()
//                                     : _robotInterface->getParams())
//      .collection.lookup(cellName);
  ControlParameterValueKind kind = parameter._kind;
  ControlParameterValue oldValue = parameter.get(kind);

  bool success = true;

  try {
    parameter.setFromString(ui->userControlTable->item(row, 1)->text().toStdString());
  } catch (std::exception& e) {
    success = false;
  }

  if (!success) {
    printf("[ERROR] invalid data, restoring old data!\n");
    parameter.set(oldValue, kind);

    assert(!_ignoreTableCallbacks);

    _ignoreTableCallbacks = true;
    ui->userControlTable->item(row, 1)->setText(
        QString(_userParameters
                    .collection.lookup(cellName)
                    .toString()
                    .c_str()));
    _ignoreTableCallbacks = false;
  } else {
    if(_started) {
      if (_simulationMode) {
        if (_simulation && _simulation->isRobotConnected()) {
          _simulation->sendControlParameter(
              cellName, parameter.get(parameter._kind), parameter._kind, true);
        }

        _ignoreTableCallbacks = true;
        ui->userControlTable->item(row, 1)->setText(
            QString(parameter.toString().c_str()));
        _ignoreTableCallbacks = false;
      } else {
      _robotInterface->sendControlParameter(
          cellName, parameter.get(parameter._kind), parameter._kind, true);
      }
    }
  }
}

void SimControlPanel::on_loadUserButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(
      nullptr, ("Load User Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  _userParameters.lockMutex();
  _userParameters.collection.deleteAll();
  _userParameters.defineAndInitializeFromYamlFile(
      fileName.toStdString());
  loadUserParameters(_userParameters);
  _userParameters.unlockMutex();
  _loadedUserSettings = true;

  if(_started) {
    if (_simulationMode) {
      if (_simulation && _simulation->isRobotConnected()) {
        for (auto& kv : _userParameters.collection._map) {
          _simulation->sendControlParameter(
              kv.first, kv.second->get(kv.second->_kind), kv.second->_kind, true);
        }
      }
    } else {
    for (auto& kv : _userParameters.collection._map) {
      _robotInterface->sendControlParameter(
          kv.first, kv.second->get(kv.second->_kind), kv.second->_kind, true);
    }

    }
  }
}

void SimControlPanel::on_saveUserButton_clicked() {
  QString fileName = QFileDialog::getSaveFileName(
      nullptr, ("Save User Table Values"), "../config", "All Files (*)");
  if (fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }
  _userParameters.writeToYamlFile(fileName.toStdString());
}

void SimControlPanel::on_goHomeButton_clicked() {
  printf("go home\n");
  FBModelState<double> homeState;
  homeState.bodyOrientation << 1, 0, 0, 0;
  homeState.bodyPosition = Vec3<double>(0, 0, 0.4);
  homeState.bodyVelocity = SVec<double>::Zero();
  homeState.q = DVec<double>(12);
  homeState.q << -0.05, -0.8, 1.7, 0.05, -0.8, 1.7, -0.05, -0.8, 1.7, 0.05, -0.8, 1.7;
  homeState.qd = homeState.q;

  _simulation->setRobotState(homeState);
}

void SimControlPanel::on_kickButton_clicked() {
  // velocity of the floating base:
  SVec<double> kickVelocity;
  kickVelocity << ui->kickAngularX->text().toDouble(),
      ui->kickAngularY->text().toDouble(), ui->kickAngularZ->text().toDouble(),
      ui->kickLinearX->text().toDouble(), ui->kickLinearY->text().toDouble(),
      ui->kickLinearZ->text().toDouble();

  FBModelState<double> state = _simulation->getRobotState();
  state.bodyVelocity += kickVelocity;
  _simulation->setRobotState(state);
}
