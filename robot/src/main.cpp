/*!
 * @file main.cpp
 * @brief Main Function for the robot program
 *
 * The main function parses command line arguments and starts the appropriate driver.
 */


#include <iostream>
#include <main.h>
#include <cassert>
#include <HardwareBridge.h>

#include "SimulationBridge.h"

MasterConfig gMasterConfig;

/*!
 * Print a message describing the command line flags for the robot program
 */
void printUsage() {
  printf("Usage: robot [robot-id] [sim-or-robot]\n\twhere robot-id:     3 for cheetah 3, m for mini-cheetah\n\t      sim-or-robot: s for sim, r for robot\n");
}

int main(int argc, char** argv) {

  if(argc != 3) {
    printUsage();
    return EXIT_FAILURE;
  }

  if(argv[1][0] == '3') {
    gMasterConfig._robot = RobotType::CHEETAH_3;
  } else if(argv[1][0] == 'm') {
    gMasterConfig._robot = RobotType::MINI_CHEETAH;
  } else {
    printUsage();
    return EXIT_FAILURE;
  }

  if(argv[2][0] == 's') {
    gMasterConfig.simulated = true;
  } else if(argv[2][0] == 'r') {
    gMasterConfig.simulated = false;
  } else {
    printUsage();
    return EXIT_FAILURE;
  }

  printf("[Quadruped] Cheetah Software\n");
  printf("        Quadruped:  %s\n", gMasterConfig._robot == RobotType::MINI_CHEETAH ? "Mini Cheetah" : "Cheetah 3");
  printf("        Driver: %s\n", gMasterConfig.simulated ? "Development Simulation Driver" : "Quadruped Driver");

  // dispatch the appropriate driver
  if(gMasterConfig.simulated) {
    if(gMasterConfig._robot == RobotType::MINI_CHEETAH) {
      SimulationBridge simulationBridge(gMasterConfig._robot);
      simulationBridge.run();
      printf("[Quadruped] SimDriver run() has finished!\n");
    } else if (gMasterConfig._robot == RobotType::CHEETAH_3) {
      SimulationBridge simulationBridge(gMasterConfig._robot);
      simulationBridge.run();
    } else {
      printf("[ERROR] unknown robot\n");
      assert(false);
    }
  } else {
    if(gMasterConfig._robot == RobotType::MINI_CHEETAH) {
      MiniCheetahHardwareBridge hw;
      hw.run();
      printf("[Quadruped] SimDriver run() has finished!\n");
    } else if (gMasterConfig._robot == RobotType::CHEETAH_3) {
      printf("[ERROR] can't do cheetah 3 hardware\n");
      assert(false);
    } else {
      printf("[ERROR] unknown robot\n");
      assert(false);
    }
  }

  return 0;
}