/*! @file main.cpp
 *  @brief Main function for simulator
 */


#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/Quadruped.h"
#include "Utilities/utilities.h"
#include "Dynamics/DynamicsSimulator.h"
#include "Dynamics/Cheetah3.h"
#include "Collision/CollisionPlane.h"
#include "Graphics3D.h"
#include "DrawList.h"
#include "Dynamics/MiniCheetah.h"
#include "Simulation.h"
#include "SimControlPanel.h"

#include <QApplication>
#include <QSurfaceFormat>

#include <unistd.h>
#include <thread>
#include <stdio.h>


/*!
 * Setup QT and run a simulation
 */
int main(int argc, char *argv[]) {

  // set up Qt
  QApplication a(argc, argv);

  // open simulator UI
  SimControlPanel panel;
  panel.show();

  // run the Qt program
  a.exec();

  return 0;
}


