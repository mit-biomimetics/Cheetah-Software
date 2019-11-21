/*============================= FSM State =============================*/
#ifndef TRANSITIONDATA_H
#define TRANSITIONDATA_H

/**
 * Struct of relevant data that can be used during transition to pass
 * data between states.
 */
template <typename T>
struct TransitionData {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  TransitionData() { zero(); }

  // Zero out all of the data
  void zero() {
    // Flag to mark when transition is done
    done = false;

    // Timing parameters
    t0 = 0.0;         // time that transition started
    tCurrent = 0.0;   // current time since transition started
    tDuration = 0.0;  // overall transition duration

    // Robot state at the beginning of transition
    comState0 = Vec12<T>::Zero();  // center of mass state
    qJoints0 = Vec12<T>::Zero();   // joint positions
    pFoot0 = Mat34<T>::Zero();     // foot positions

    // Current robot state
    comState = Vec12<T>::Zero();  // center of mass state
    qJoints = Vec12<T>::Zero();   // joint positions
    pFoot = Mat34<T>::Zero();     // foot positions
  }

  // Flag to mark when transition is done
  bool done = false;

  // Timing parameters
  T t0;         // time that transition started
  T tCurrent;   // current time since transition started
  T tDuration;  // overall transition duration

  // Robot state at the beginning of transition
  Vec12<T> comState0;  // center of mass state
  Vec12<T> qJoints0;   // joint positions
  Mat34<T> pFoot0;     // foot positions

  // Current robot state
  Vec12<T> comState;  // center of mass state
  Vec12<T> qJoints;   // joint positions
  Mat34<T> pFoot;     // foot positions
};

template struct TransitionData<double>;
template struct TransitionData<float>;

#endif  // CONTROLFSM_H