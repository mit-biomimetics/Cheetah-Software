#ifndef JOYSTICK_TEST_H
#define JOYSTICK_TEST_H

#include <QDialog>
#include <SimUtilities/GamepadCommand.h>
#include "GameController.h"

namespace Ui {
  class JoystickTestWindow;
}

class JoystickTestWindow : public QDialog {
Q_OBJECT

public:
  explicit JoystickTestWindow(GameController &gamepad, QWidget *parent = nullptr);

  ~JoystickTestWindow();

private:
  Ui::JoystickTestWindow *ui;
  GameController& _gamepad;
  GamepadCommand _command;

private slots:
  void update();
};

#endif // JOYSTICK_TEST_H
