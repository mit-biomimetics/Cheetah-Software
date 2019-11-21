#include "JPos_Controller.hpp"

void JPos_Controller::runController(){

  Mat3<float> kpMat;
  Mat3<float> kdMat;
  //kpMat << 20, 0, 0, 0, 20, 0, 0, 0, 20;
  //kdMat << 2.1, 0, 0, 0, 2.1, 0, 0, 0, 2.1;
  kpMat << userParameters.kp, 0, 0, 0,  userParameters.kp, 0, 0, 0,  userParameters.kp;
  kdMat <<  userParameters.kd, 0, 0, 0, userParameters.kd, 0, 0, 0, userParameters.kd;

  static int iter(0);
  ++iter;

  if(iter < 10){
    for(int leg(0); leg<4; ++leg){
      for(int jidx(0); jidx<3; ++jidx){
        _jpos_ini[3*leg+jidx] = _legController->datas[leg].q[jidx];
      }
    }
  }

  _legController->_maxTorque = 150;
  _legController->_legsEnabled = true;

  if(userParameters.calibrate > 0.4) {
    _legController->_calibrateEncoders = userParameters.calibrate;
  } else {
    if(userParameters.zero > 0.5) {
      _legController->_zeroEncoders = true;
    } else {
      _legController->_zeroEncoders = false;

      for(int leg(0); leg<4; ++leg){
        for(int jidx(0); jidx<3; ++jidx){
          float pos = std::sin(.001f * iter);
          _legController->commands[leg].qDes[jidx] = pos;
          _legController->commands[leg].qdDes[jidx] = 0.;
          _legController->commands[leg].tauFeedForward[jidx] = userParameters.tau_ff;
        }
        _legController->commands[leg].kpJoint = kpMat;
        _legController->commands[leg].kdJoint = kdMat;
      }
    }
  }



  //if(iter%200 ==0){
    //printf("value 1, 2: %f, %f\n", userParameters.testValue, userParameters.testValue2);
  //}


}
