//#include "HardwareInterface.h"
//#include "main.h"
//#include <sys/timerfd.h>
//#include <thread>
//#include <pthread.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
//#include <iostream>
//#include <stdio.h>
//#include <sched.h>
//#include <stdlib.h>
//#include <sys/time.h>
//#include <sys/mman.h>
//#include <unistd.h>
//#include <signal.h>
//#include <string.h>
//#include <pthread.h>
//#include <time.h>
//#include <lcm/lcm.h>
//#include <thread>
//
// static constexpr int taskPriority = 49;
//
///*!
// * Initialize the robot hardware
// */
// HardwareInterface::HardwareInterface() {
//  prefaultStack();
//  setupRT();
//
//}
//
///*!
// * Make sure there are no page faults when accessing the stack
// */
// void HardwareInterface::prefaultStack() {
//  constexpr int stackSize = 24*1024;
//  volatile unsigned char stack[stackSize];
//  stack[stackSize/2] = 2;
//  memset((void*)stack,0,stackSize);
//  if(mlockall(MCL_CURRENT|MCL_FUTURE) == -1) {
//    printf("[HardwareInterface] failed to mlockall.  This is likely because
//    cheetah wasn't run as root, and is okay for simulation.\n");
//  }
//}
//
// void HardwareInterface::setupRT() {
//    sched_param params;
//    params.sched_priority = taskPriority;
//    sched_setscheduler(0, SCHED_FIFO, &params);
//}