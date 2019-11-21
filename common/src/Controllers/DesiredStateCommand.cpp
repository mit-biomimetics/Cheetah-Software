/*========================= Gamepad Control ==========================*/
/**
 *
 */

#include "Controllers/DesiredStateCommand.h"
/*=========================== Gait Data ===============================*/
/**
 *
 */
template <typename T>
void DesiredStateData<T>::zero() {
  // Overall desired state
  stateDes = Vec12<T>::Zero();
  stateTrajDes = Eigen::Matrix<T, 12, 10>::Zero();
}

template struct DesiredStateData<double>;
template struct DesiredStateData<float>;

/**
 *
 */
template <typename T>
void DesiredStateCommand<T>::convertToStateCommands() {
  data.zero();
  Vec2<float> joystickLeft, joystickRight;

  //T height_cmd(0.3);

  // THIS SHOULD be DISABLE Soon
  if(parameters->use_rc) {
    if(rcCommand->mode == RC_mode::QP_STAND){ // Stand
      joystickLeft[0] = 0.; // Y
      joystickLeft[1] = 0.;
      joystickRight[0] = rcCommand->rpy_des[2]; // Yaw
      //height_cmd = rcCommand->height_variation;

    }else if(rcCommand->mode == RC_mode::LOCOMOTION ||
        rcCommand->mode == RC_mode::VISION){ // Walking
      joystickLeft[0] = rcCommand->v_des[1]; // Y
      joystickLeft[1] = rcCommand->v_des[0]; // X
      joystickRight[0] = rcCommand->omega_des[2]; // Yaw
      joystickRight[1] = rcCommand->omega_des[1]; // Pitch
      //height_cmd = rcCommand->height_variation;

    }else if(rcCommand->mode == RC_mode::TWO_LEG_STANCE){ // Two Contact Stand
      //joystickLeft[0] = rcCommand->p_des[1]; // Y
      joystickLeft[1] = rcCommand->p_des[0]; // X
      joystickRight[0] = rcCommand->rpy_des[2]; // Yaw
      joystickRight[1] = rcCommand->rpy_des[1]; // Pitch
      joystickLeft[0] = rcCommand->rpy_des[0]; // Roll

    }else{
      joystickLeft.setZero();
      joystickRight.setZero();
    }
  } else { // No Remote Controller
    joystickLeft = gamepadCommand->leftStickAnalog;
    joystickRight = gamepadCommand->rightStickAnalog;
    trigger_pressed = gamepadCommand->a;
  }
  // Warning!!!!
  // Recommend not to use stateDes
  // We are going to remove it soon

  joystickLeft[0] *= -1.f;
  joystickRight[0] *= -1.f;

  leftAnalogStick = leftAnalogStick * (T(1) - filter) + joystickLeft * filter;
  rightAnalogStick = rightAnalogStick * (T(1) - filter) + joystickRight * filter;

  // Desired states from the controller
  data.stateDes(6) = deadband(leftAnalogStick[1], minVelX, maxVelX);  // forward linear velocity
  data.stateDes(7) = deadband(leftAnalogStick[0], minVelY, maxVelY);  // lateral linear velocity
  data.stateDes(8) = 0.0;  // vertical linear velocity
  data.stateDes(0) = dt * data.stateDes(6);  // X position
  data.stateDes(1) = dt * data.stateDes(7);  // Y position
  data.stateDes(2) = 0.26;  // Z position height
  data.stateDes(9) = 0.0;  // Roll rate
  data.stateDes(10) = 0.0;  // Pitch rate
  data.stateDes(11) = deadband(rightAnalogStick[0], minTurnRate, maxTurnRate);  // Yaw turn rate
  data.stateDes(3) = 0.0; // Roll
  data.stateDes(4) = deadband(rightAnalogStick[1], minPitch, maxPitch);  // Pitch
  data.stateDes(5) = dt * data.stateDes(11);  // Yaw
}

template <typename T>
void DesiredStateCommand<T>::setCommandLimits(T minVelX_in, T maxVelX_in,
    T minVelY_in, T maxVelY_in, T minTurnRate_in, T maxTurnRate_in) {
  minVelX = minVelX_in;
  maxVelX = maxVelX_in;
  minVelY = minVelY_in;
  maxVelY = maxVelY_in;
  minTurnRate = minTurnRate_in;
  maxTurnRate = maxTurnRate_in;
}

/**
 *
 */
template <typename T>
void DesiredStateCommand<T>::desiredStateTrajectory(int N, Vec10<T> dtVec) {
  A = Mat12<T>::Zero();
  A(0, 0) = 1;
  A(1, 1) = 1;
  A(2, 2) = 1;
  A(3, 3) = 1;
  A(4, 4) = 1;
  A(5, 5) = 1;
  A(6, 6) = 1;
  A(7, 7) = 1;
  A(8, 8) = 1;
  A(9, 9) = 1;
  A(10, 10) = 1;
  A(11, 11) = 1;
  data.stateTrajDes.col(0) = data.stateDes;

  for (int k = 1; k < N; k++) {
    A(0, 6) = dtVec(k - 1);
    A(1, 7) = dtVec(k - 1);
    A(2, 8) = dtVec(k - 1);
    A(3, 9) = dtVec(k - 1);
    A(4, 10) = dtVec(k - 1);
    A(5, 11) = dtVec(k - 1);
    data.stateTrajDes.col(k) = A * data.stateTrajDes.col(k - 1);
    for (int i = 0; i < 12; i++) {
      // std::cout << data.stateTrajDes(i, k) << " ";
    }
    // std::cout << std::endl;
  }
  // std::cout << std::endl;
}

/**
 *
 */
template <typename T>
float DesiredStateCommand<T>::deadband(float command, T minVal, T maxVal) {
  if (command < deadbandRegion && command > -deadbandRegion) {
    return 0.0;
  } else {
    return (command / (2)) * (maxVal - minVal);
  }
}

/**
 *
 */
template <typename T>
void DesiredStateCommand<T>::printRawInfo() {
  // Increment printing iteration
  printIter++;

  // Print at requested frequency
  if (printIter == printNum) {
    std::cout << "[DESIRED STATE COMMAND] Printing Raw Gamepad Info...\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "Button Start: " << gamepadCommand->start
              << " | Back: " << gamepadCommand->back << "\n";
    std::cout << "Button A: " << gamepadCommand->a
              << " | B: " << gamepadCommand->b << " | X: " << gamepadCommand->x
              << " | Y: " << gamepadCommand->y << "\n";
    std::cout << "Left Stick Button: " << gamepadCommand->leftStickButton
              << " | X: " << gamepadCommand->leftStickAnalog[0]
              << " | Y: " << gamepadCommand->leftStickAnalog[1] << "\n";
    std::cout << "Right Analog Button: " << gamepadCommand->rightStickButton
              << " | X: " << gamepadCommand->rightStickAnalog[0]
              << " | Y: " << gamepadCommand->rightStickAnalog[1] << "\n";
    std::cout << "Left Bumper: " << gamepadCommand->leftBumper
              << " | Trigger Switch: " << gamepadCommand->leftTriggerButton
              << " | Trigger Value: " << gamepadCommand->leftTriggerAnalog
              << "\n";
    std::cout << "Right Bumper: " << gamepadCommand->rightBumper
              << " | Trigger Switch: " << gamepadCommand->rightTriggerButton
              << " | Trigger Value: " << gamepadCommand->rightTriggerAnalog
              << "\n\n";
    std::cout << std::endl;

    // Reset iteration counter
    printIter = 0;
  }
}

/**
 *
 */
template <typename T>
void DesiredStateCommand<T>::printStateCommandInfo() {
  // Increment printing iteration
  printIter++;

  // Print at requested frequency
  if (printIter == printNum) {
    std::cout << "[DESIRED STATE COMMAND] Printing State Command Info...\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "Position X: " << data.stateDes(0)
              << " | Y: " << data.stateDes(1) << " | Z: " << data.stateDes(2)
              << "\n";
    std::cout << "Orientation Roll: " << data.stateDes(3)
              << " | Pitch: " << data.stateDes(4)
              << " | Yaw: " << data.stateDes(5) << "\n";
    std::cout << "Velocity X: " << data.stateDes(6)
              << " | Y: " << data.stateDes(7) << " | Z: " << data.stateDes(8)
              << "\n";
    std::cout << "Angular Velocity X: " << data.stateDes(9)
              << " | Y: " << data.stateDes(10) << " | Z: " << data.stateDes(11)
              << "\n";
    std::cout << std::endl;
    std::cout << std::endl;

    // Reset iteration counter
    printIter = 0;
  }
}

template class DesiredStateCommand<double>;
template class DesiredStateCommand<float>;
