#ifndef PROJECT_TYPES_H
#define PROJECT_TYPES_H

#include "cppTypes.h"

struct MasterConfig {
  RobotType _robot;
  bool simulated = false;
  bool load_from_file = false;
};

#endif  // PROJECT_TYPES_H
