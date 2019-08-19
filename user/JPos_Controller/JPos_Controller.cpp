#include "JPos_Controller.hpp"

void JPos_Controller::runController(){

  Mat3<float> kpMat;
  Mat3<float> kdMat;
  kpMat << 20, 0, 0, 0, 20, 0, 0, 0, 20;
  kdMat << 2.1, 0, 0, 0, 2.1, 0, 0, 0, 2.1;
  //kpMat << 1., 0, 0, 0, 1., 0, 0, 0, 1.;
  //kdMat << 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01;

  static int iter(0);
  ++iter;

  if(iter < 10){
    for(int leg(0); leg<4; ++leg){
      for(int jidx(0); jidx<3; ++jidx){
        _jpos_ini[3*leg+jidx] = _legController->datas[leg].q[jidx];
      }
    }
  }

  for(int leg(0); leg<4; ++leg){
    for(int jidx(0); jidx<3; ++jidx){
      _legController->commands[leg].qDes[jidx] = _jpos_ini[3*leg + jidx];
      _legController->commands[leg].qdDes[jidx] = 0.;
      _legController->commands[leg].tauFeedForward[jidx] = 0.;
    }
    _legController->commands[leg].kpJoint = kpMat;
    _legController->commands[leg].kdJoint = kdMat;

    _legController->commands[leg].qDes[2] += 0.6 * sin(iter * 0.002);
  }
}
