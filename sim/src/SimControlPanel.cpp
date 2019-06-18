#include "SimControlPanel.h"
#include <QMessageBox>
#include <QFileDialog>
#include <ControlParameters/ControlParameters.h>
#include "ui_SimControlPanel.h"

/*!
 * Display an error messagebox with the given text
 */
static void createErrorMessage(const std::string& text) {
  QMessageBox mb;
  mb.setText(QString(text.c_str()));
  mb.exec();
}


SimControlPanel::SimControlPanel(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SimControlPanel),
    _terrainFileName(getConfigDirectoryPath() + DEFAULT_TERRAIN_FILE){
  ui->setupUi(this);
  updateUiEnable();
  updateTerrainLabel();

  printf("[SimControlPanel] Init simulator parameters...\n");
  _parameters.initializeFromYamlFile(getConfigDirectoryPath() + SIMULATOR_DEFAULT_PARAMETERS);
  if(!_parameters.isFullyInitialized()) {
    printf("[ERROR] Simulator parameters are not fully initialized.  You forgot: \n%s\n", _parameters.generateUnitializedList().c_str());
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

void SimControlPanel::updateTerrainLabel() {
  ui->terrainFileLabel->setText(QString(_terrainFileName.c_str()));
}

void SimControlPanel::on_startButton_clicked() {
  RobotType robotType;

  if(ui->cheetah3Button->isChecked()) {
    robotType = RobotType::CHEETAH_3;
  } else if(ui->miniCheetahButton->isChecked()) {
    robotType = RobotType::MINI_CHEETAH;
  } else {
    createErrorMessage("Error: you must select a robot");
    return;
  }

  if(!ui->simulatorButton->isChecked() && !ui->robotButton->isChecked()) {
    createErrorMessage("Error: you must select either robot or simulation mode");
    return;
  }

  _simulationMode = ui->simulatorButton->isChecked();

  printf("[SimControlPanel] Initialize Graphics...\n");
  _graphicsWindow = new Graphics3D();
  _graphicsWindow->show();
  _graphicsWindow->resize(1280, 720);

  if(_simulationMode) {
    printf("[SimControlPanel] Initialize simulator...\n");
    _simulation = new Simulation(robotType, _graphicsWindow, _parameters);
    loadSimulationParameters(_simulation->getSimParams());
    loadRobotParameters(_simulation->getRobotParams());

    printf("[SimControlParameter] Load terrain...\n");
    _simulation->loadTerrainFile(_terrainFileName);
    _simThread = std::thread([this](){_simulation->runAtSpeed();});

    _graphicsWindow->setAnimating(true);
  } else {
    printf("[SimControlPanel] Init Robot Interface...\n");
    _interfaceTaskManager = new PeriodicTaskManager;
    _robotInterface = new RobotInterface(robotType, _graphicsWindow, _interfaceTaskManager);
    loadRobotParameters(_robotInterface->getParams());
    _robotInterface->startInterface();
    _graphicsWindow->setAnimating(true);
  }

  _started = true;
  updateUiEnable();
}



void SimControlPanel::on_stopButton_clicked() {
  if(_simulation) {
    _simulation->stop();
    _simThread.join();
  } else {
    _robotInterface->stopInterface();
  }

  if(_graphicsWindow) {
    _graphicsWindow->setAnimating(false);
    _graphicsWindow->hide();
  }

  printf("calling destructors\n");
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

static void updateQtableWithParameters(ControlParameters& params, QTableWidget& table) {
  table.setRowCount((s32)params.collection._map.size());

  s32 i = 0;
  for(auto& kv : params.collection._map) {
    for(s32 col = 0; col < 2; col++) {
      QTableWidgetItem* cell = table.item(i,col);
      if(!cell) {
        cell = new QTableWidgetItem;
        table.setItem(i,col,cell);
      }
    }

    table.item(i,0)->setText(QString(kv.first.c_str()));
    table.item(i,1)->setText(QString(kv.second->toString().c_str()));
    i++;
  }
}

/*!
 * Reload all values in the Simulation Parameter Table
 */
void SimControlPanel::loadSimulationParameters(SimulatorControlParameters &params) {
  _ignoreTableCallbacks = true;
  updateQtableWithParameters(params, *ui->simulatorTable);
  _ignoreTableCallbacks = false;
}

void SimControlPanel::loadRobotParameters(RobotControlParameters &params) {
  _ignoreTableCallbacks = true;
  updateQtableWithParameters(params, *ui->robotTable);
  _ignoreTableCallbacks = false;
}

void SimControlPanel::on_joystickButton_clicked() {
  _graphicsWindow->resetGameController();
}

void SimControlPanel::on_driverButton_clicked() {

}

void SimControlPanel::on_simulatorTable_cellChanged(int row, int column) {
  if(_ignoreTableCallbacks) return;

  if(column != 1) {
    return;
  }

  auto cell = ui->simulatorTable->item(row, 0);
  std::string cellName = cell->text().toStdString();

  if(cellName == "") {
    return;
  }



  auto& parameter = _parameters.collection.lookup(cellName);
  ControlParameterValueKind kind = parameter._kind;
  ControlParameterValue oldValue = parameter.get(kind);

  bool success = true;

  try {
    parameter.setFromString(ui->simulatorTable->item(row, 1)->text().toStdString());
  } catch (std::exception& e) {
    success = false;
  }

  if(!success) {
    printf("[ERROR] invalid data, restoring old data!\n");
    parameter.set(oldValue, kind);

    assert(!_ignoreTableCallbacks);

    _ignoreTableCallbacks = true;
    ui->simulatorTable->item(row, 1)->setText(QString(_parameters.collection.lookup(cellName).toString().c_str()));
    _ignoreTableCallbacks = false;
  } else {
    // this update "rewrites" the value.  If it's an integer, it kills any decimal.  If it's a float, it puts in
    // scientific notation if needed.
    _ignoreTableCallbacks = true;
    ui->simulatorTable->item(row, 1)->setText(QString(parameter.toString().c_str()));
    _ignoreTableCallbacks = false;
  }


}

void SimControlPanel::on_saveSimulatorButton_clicked() {
  QString fileName = QFileDialog::getSaveFileName(nullptr, ("Save Simulator Table Values"), "../config", "All Files (*)");
  if(fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }


  _parameters.lockMutex();
  _parameters.writeToYamlFile(fileName.toStdString());
  _parameters.unlockMutex();
}

void SimControlPanel::on_loadSimulatorButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(nullptr, ("Load Simulator Table Values"), "../config", "All Files (*)");
  if(fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }
  ;
  _parameters.collection.clearAllSet();
  _parameters.initializeFromYamlFile(fileName.toStdString());
  if(!_parameters.collection.checkIfAllSet()) {
    printf("new settings file %s doesn't contain the following simulator parameters:\n%s\n",
            fileName.toStdString().c_str(), _parameters.generateUnitializedList().c_str());
    throw std::runtime_error("bad new settings file");
  }
  loadSimulationParameters(_parameters);
  _parameters.unlockMutex();
}

void SimControlPanel::on_robotTable_cellChanged(int row, int column) {
  if(_ignoreTableCallbacks) return;
  if(column != 1) {
    return;
  }

  auto cell = ui->robotTable->item(row, 0);
  std::string cellName = cell->text().toStdString();

  if(cellName == "") {
    return;
  }

  auto& parameter = (_simulationMode ? _simulation->getRobotParams() : _robotInterface->getParams()).collection.lookup(cellName);
  ControlParameterValueKind kind = parameter._kind;
  ControlParameterValue oldValue = parameter.get(kind);

  bool success = true;

  try {
    parameter.setFromString(ui->robotTable->item(row, 1)->text().toStdString());
  } catch (std::exception& e) {
    success = false;
  }

  if(!success) {
    printf("[ERROR] invalid data, restoring old data!\n");
    parameter.set(oldValue, kind);

    assert(!_ignoreTableCallbacks);

    _ignoreTableCallbacks = true;
    ui->robotTable->item(row, 1)->setText(QString(_simulation->getRobotParams().collection.lookup(cellName).toString().c_str()));
    _ignoreTableCallbacks = false;
  } else {
    if(_simulationMode) {
      if(_simulation->isRobotConnected()) {
        _simulation->sendControlParameter(cellName, parameter.get(parameter._kind), parameter._kind);
      }

      _ignoreTableCallbacks = true;
      ui->robotTable->item(row, 1)->setText(QString(parameter.toString().c_str()));
      _ignoreTableCallbacks = false;
    } else {
      _robotInterface->sendControlParameter(cellName, parameter.get(parameter._kind), parameter._kind);
    }
  }




}

void SimControlPanel::on_saveRobotButton_clicked() {
  printf("save callback\n");
  QString fileName = QFileDialog::getSaveFileName(nullptr, ("Save Quadruped Table Values"), "../config", "All Files (*)");
  printf("2\n");
  if(fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }
  printf("3\n");
  _simulation->getRobotParams().writeToYamlFile(fileName.toStdString());
}

void SimControlPanel::on_loadRobotButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(nullptr, ("Load Quadruped Table Values"), "../config", "All Files (*)");
  if(fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  if(_simulationMode) {
    _simulation->getRobotParams().lockMutex();
    _simulation->getRobotParams().collection.clearAllSet();
    _simulation->getRobotParams().initializeFromYamlFile(fileName.toStdString());
    if(!_simulation->getRobotParams().collection.checkIfAllSet()) {
      printf("new settings file %s doesn't contain the following robot parameters:\n%s\n",
             fileName.toStdString().c_str(), _simulation->getRobotParams().generateUnitializedList().c_str());
      throw std::runtime_error("bad new settings file");
    }
    loadRobotParameters(_simulation->getRobotParams());

    if(_simulation->isRobotConnected()) {
      for(auto& kv : _simulation->getRobotParams().collection._map) {
        _simulation->sendControlParameter(kv.first, kv.second->get(kv.second->_kind), kv.second->_kind);
      }
    }
    _simulation->getSimParams().unlockMutex();
  } else {
    _robotInterface->getParams().lockMutex();
    _robotInterface->getParams().collection.clearAllSet();
    _robotInterface->getParams().initializeFromYamlFile(fileName.toStdString());
    if(!_robotInterface->getParams().collection.checkIfAllSet()) {
      printf("new settings file %s doesn't contain the following robot parameters:\n%s\n",
             fileName.toStdString().c_str(), _robotInterface->getParams().generateUnitializedList().c_str());
      throw std::runtime_error("bad new settings file");
    }
    loadRobotParameters(_robotInterface->getParams());

    for(auto& kv : _robotInterface->getParams().collection._map) {
      _robotInterface->sendControlParameter(kv.first, kv.second->get(kv.second->_kind), kv.second->_kind);
    }

    _robotInterface->getParams().unlockMutex();
  }
}

void SimControlPanel::on_setTerrainButton_clicked() {
  QString fileName = QFileDialog::getOpenFileName(nullptr, ("Load Terrain Definition"), "../config", "All Files (*)");
  if(fileName == nullptr || fileName == "") {
    createErrorMessage("File name is invalid");
    return;
  }

  _terrainFileName = fileName.toStdString();
  updateTerrainLabel();
}

void SimControlPanel::on_favoritesTable_cellChanged(int row, int column) {
  (void) row;
  (void) column;
}

void SimControlPanel::on_loadFavoriteButton_clicked() {

}


void SimControlPanel::on_goHomeButton_clicked() {
  printf("go home\n");
  FBModelState<double> homeState;
  homeState.bodyOrientation << 1, 0, 0, 0;
  homeState.bodyPosition = Vec3<double>(0,0,0.5);
  homeState.bodyVelocity = SVec<double>::Zero();
  homeState.q = DVec<double>(12);
  homeState.q << 0,0,0,0,0,0,0,0,0,0,0,0;
  homeState.qd = homeState.q;

  _simulation->setRobotState(homeState);
}

void SimControlPanel::on_kickButton_clicked() {
  // velocity of the floating base:
  SVec<double> kickVelocity;
  kickVelocity <<
  ui->kickAngularX->text().toDouble(),
  ui->kickAngularY->text().toDouble(),
  ui->kickAngularZ->text().toDouble(),
  ui->kickLinearX->text().toDouble(),
  ui->kickLinearY->text().toDouble(),
  ui->kickLinearZ->text().toDouble();

  FBModelState<double> state = _simulation->getRobotState();
  state.bodyVelocity += kickVelocity;
  _simulation->setRobotState(state);
}
