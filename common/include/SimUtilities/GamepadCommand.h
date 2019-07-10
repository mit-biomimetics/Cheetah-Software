/*! @file GamepadCommand.h
 *  @brief The GamepadCommand type containing joystick information
 */

#ifndef PROJECT_GAMEPADCOMMAND_H
#define PROJECT_GAMEPADCOMMAND_H

#include "Utilities/utilities.h"
#include "cppTypes.h"
#include "gamepad_lcmt.hpp"

struct GamepadCommand {
  GamepadCommand() { zero(); }

  bool leftBumper, rightBumper, leftTriggerButton, rightTriggerButton, back,
      start, a, b, x, y, leftStickButton, rightStickButton, logitechButton;

  Vec2<float> leftStickAnalog, rightStickAnalog;
  float leftTriggerAnalog, rightTriggerAnalog;

  void zero() {
    leftBumper = false;
    rightBumper = false;
    leftTriggerButton = false;
    rightTriggerButton = false;
    back = false;
    start = false;
    a = false;
    b = false;
    x = false;
    y = false;
    leftStickButton = false;
    rightStickButton = false;

    leftTriggerAnalog = 0;
    rightTriggerAnalog = 0;
    leftStickAnalog = Vec2<float>::Zero();
    rightStickAnalog = Vec2<float>::Zero();
  }

  void set(const gamepad_lcmt* lcmt) {
    leftBumper = lcmt->leftBumper;
    rightBumper = lcmt->rightBumper;
    leftTriggerButton = lcmt->leftTriggerButton;
    rightTriggerButton = lcmt->rightTriggerButton;
    back = lcmt->back;
    start = lcmt->start;
    a = lcmt->a;
    x = lcmt->x;
    b = lcmt->b;
    y = lcmt->y;
    leftStickButton = lcmt->leftStickButton;
    rightStickButton = lcmt->rightStickButton;
    leftTriggerAnalog = lcmt->leftTriggerAnalog;
    rightTriggerAnalog = lcmt->rightTriggerAnalog;
    for (int i = 0; i < 2; i++) {
      leftStickAnalog[i] = lcmt->leftStickAnalog[i];
      rightStickAnalog[i] = lcmt->rightStickAnalog[i];
    }
  }

  void get(gamepad_lcmt* lcmt) {
    lcmt->leftBumper = leftBumper;
    lcmt->rightBumper = rightBumper;
    lcmt->leftTriggerButton = leftTriggerButton;
    lcmt->rightTriggerButton = rightTriggerButton;
    lcmt->back = back;
    lcmt->start = start;
    lcmt->a = a;
    lcmt->x = x;
    lcmt->b = b;
    lcmt->y = y;
    lcmt->leftStickButton = leftStickButton;
    lcmt->rightStickButton = rightStickButton;
    lcmt->leftTriggerAnalog = leftTriggerAnalog;
    lcmt->rightTriggerAnalog = rightTriggerAnalog;
    for (int i = 0; i < 2; i++) {
      lcmt->leftStickAnalog[i] = leftStickAnalog[i];
      lcmt->rightStickAnalog[i] = rightStickAnalog[i];
    }
  }

  /*!
   * The Logitech F310's seem to do a bad job of returning to zero exactly, so a
   * deadband around zero is useful when integrating joystick commands
   */
  void applyDeadband(float f) {
    eigenDeadband(leftStickAnalog, f);
    eigenDeadband(rightStickAnalog, f);
    leftTriggerAnalog = deadband(leftTriggerAnalog, f);
    rightTriggerAnalog = deadband(rightTriggerAnalog, f);
  }

  std::string toString() {
    std::string result =
        "Result:\nleftBumper: " + boolToString(leftBumper) + "\n" +
        "rightBumper: " + boolToString(rightBumper) + "\n" +
        "leftTriggerButton: " + boolToString(leftTriggerButton) + "\n" +
        "rightTriggerButton: " + boolToString(rightTriggerButton) + "\n" +
        "back: " + boolToString(back) + "\n" + "start: " + boolToString(start) +
        "\n" + "a: " + boolToString(a) + "\n" + "b: " + boolToString(b) + "\n" +
        "x: " + boolToString(x) + "\n" + "y: " + boolToString(y) + "\n" +
        "leftStickButton: " + boolToString(leftStickButton) + "\n" +
        "rightStickButton: " + boolToString(rightStickButton) + "\n" +
        "leftTriggerAnalog: " + std::to_string(leftTriggerAnalog) + "\n" +
        "rightTriggerAnalog: " + std::to_string(rightTriggerAnalog) + "\n" +
        "leftStickAnalog: " + eigenToString(leftStickAnalog) + "\n" +
        "rightStickAnalog: " + eigenToString(rightStickAnalog) + "\n";
    return result;
  }
};

#endif  // PROJECT_DRIVERCOMMAND_H
