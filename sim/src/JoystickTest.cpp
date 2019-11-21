#include "../include/JoystickTest.h"
#include "ui_JoystickTest.h"
#include <QTimer>


JoystickTestWindow::JoystickTestWindow(GameController& gamepad, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::JoystickTestWindow),
    _gamepad(gamepad)
{
    ui->setupUi(this);

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(update()));
  timer->start(1000 / 30);
}

JoystickTestWindow::~JoystickTestWindow()
{
    delete ui;
}


void JoystickTestWindow::update() {
  _gamepad.updateGamepadCommand(_command);
  char buffer[256];

  sprintf(buffer, "Left X: %4.2f\n", _command.leftStickAnalog[0]);
  ui->leftXLabel->setText(buffer);

  sprintf(buffer, "Left Y: %4.2f\n", _command.leftStickAnalog[1]);
  ui->leftYLabel->setText(buffer);

  sprintf(buffer, "Right X: %4.2f\n", _command.rightStickAnalog[0]);
  ui->rightXLabel->setText(buffer);

  sprintf(buffer, "Right Y: %4.2f\n", _command.rightStickAnalog[1]);
  ui->rightYLabel->setText(buffer);
}
