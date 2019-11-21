/*!
 * @file main_helper.h
 * @brief Function which should be called in main to start your robot control code
 */

#ifndef ROBOT_MAIN_H
#define ROBOT_MAIN_H

#include "Types.h"
#include "RobotController.h"

extern MasterConfig gMasterConfig;
int main_helper(int argc, char** argv, RobotController* ctrl);

#endif  // ROBOT_MAIN_H
