/*! @file GameController.h
 *  @brief Code to read the Logitech F310 Game Controller
 *  Creates a DriverCommand object to be sent to the robot controller
 *  Used in the development simulator and in the robot control mode
 */

#ifndef PROJECT_GAMECONTROLLER_H
#define PROJECT_GAMECONTROLLER_H

#include "SimUtilities/GamepadCommand.h"

#include <QtCore/QObject>

class QGamepad;  // for an unknown reason, #including <QtGamepad/QGamepad> here
                 // makes compilation *very* slow

class GameController : public QObject {
  Q_OBJECT
 public:
  explicit GameController(QObject *parent = 0);
  void updateGamepadCommand(GamepadCommand &gamepadCommand);
  void findNewController();
  ~GameController();

 private:
  QGamepad *_qGamepad = nullptr;
};

#endif  // PROJECT_GAMECONTROLLER_H
