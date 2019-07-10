#ifndef GAIT_SCHEDULER_H
#define GAIT_SCHEDULER_H

#include <cppTypes.h>
#include <iostream>
#include "Dynamics/Quadruped.h"

/**
 * Enumerated gait types. Preplanned gaits are defined.
 */
enum class GaitType {
  STAND,
  STAND_CYCLE,
  STATIC_WALK,
  AMBLE,
  TROT_WALK,
  TROT,
  TROT_RUN,
  PACE,
  BOUND,
  ROTARY_GALLOP,
  TRAVERSE_GALLOP,
  PRONK,
  THREE_FOOT,
  CUSTOM,
  TRANSITION_TO_STAND
};

/**
 *
 */
template <typename T>
struct GaitData {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  GaitData() { zero(); }

  // Zero out all of the data
  void zero();

  // The current GaitType
  GaitType _currentGait;

  // Next GaitType to transition into
  GaitType _nextGait;

  // Gait name string
  std::string gaitName;

  // Gait descriptors
  T periodTimeNominal;      // overall period time to scale
  T initialPhase;           // initial phase to offset
  T switchingPhaseNominal;  // nominal phase to switch contacts

  // Enable flag for each foot
  Eigen::Vector4i gaitEnabled;  // enable gaint controlled legs

  // Time based descriptors
  Vec4<T> periodTime;           // overall foot scaled gait period time
  Vec4<T> timeStance;           // total stance time
  Vec4<T> timeSwing;            // total swing time
  Vec4<T> timeStanceRemaining;  // stance time remaining
  Vec4<T> timeSwingRemaining;   // swing time remaining

  // Phase based descriptors
  Vec4<T> switchingPhase;  // phase to switch to swing
  Vec4<T> phaseVariable;   // overall gait phase for each foot
  Vec4<T> phaseOffset;     // nominal gait phase offsets
  Vec4<T> phaseScale;      // phase scale relative to variable
  Vec4<T> phaseStance;     // stance subphase
  Vec4<T> phaseSwing;      // swing subphase

  // Scheduled contact states
  Eigen::Vector4i contactStateScheduled;  // contact state of the foot
  Eigen::Vector4i contactStatePrev;       // previous contact state of the foot
  Eigen::Vector4i touchdownScheduled;     // scheduled touchdown event flag
  Eigen::Vector4i liftoffScheduled;       // scheduled touchdown event flag

  // Position of the feet in the world frame at takeoff time
  Mat34<T> posFootLiftoffWorld;    // foot position when scheduled to lift off
  Mat34<T> posFootTouchdownWorld;  // foot position when scheduled to touchdown
};

/**
 *
 */
template <typename T>
class GaitScheduler {
 public:
  // Constructors for the GaitScheduler
  GaitScheduler();
  ~GaitScheduler(){};

  // Initialize the Gait Scheduler
  void initialize();

  // Iteration step for scheduler logic
  void step();

  // Creates a new gait from predefined library
  void createGait();

  // Prints the characteristic info and curret state
  void printGaitInfo();

  // Struct containing all of the gait relevant data
  GaitData<T> gaitData;

 private:
  // The quadruped model
  // Quadruped<T>& _quadruped;

  // Control loop timestep change
  T dt = 0.001;

  // Phase change at each step
  T dphase;

  // Choose how often to print info, every N iterations
  int printNum = 5;  // N*(0.001s) in simulation time

  // Track the number of iterations since last info print
  int printIter = 0;
};

#endif