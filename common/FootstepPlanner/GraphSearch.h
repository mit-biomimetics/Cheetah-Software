#ifndef CHEETAH_SOFTWARE_GRAPHSEARCH_H
#define CHEETAH_SOFTWARE_GRAPHSEARCH_H

#include <vector>
#include "cppTypes.h"


struct ContactState {
  union {
    bool contact[4];
    struct {
      bool fr, fl, rr, rl;
    };
  };

  ContactState(bool _fr, bool _fl, bool _rr, bool _rl) {
    fr = _fr;
    fl = _fl;
    rr = _rr;
    rl = _rl;
  }

  ContactState() { }
};

struct DefaultGaits {
  std::vector<ContactState> trotting, standing;
};

struct InputTrajectoryState {
  Vec2<float> p;
  Vec2<float> v;
  float theta;
};


struct FootplanFootState {
  Vec2<float> p;
  bool contact;
  float stateTime;
};

struct FootplanState {
  float t;
  Vec2<float> pBase;
  FootplanFootState feet[4];
};

struct FootplanStats {
  u64 nodesVisited;
  u64 maxMemory;

  FootplanStats() {
    reset();
  }

  void reset() {
    nodesVisited = 0;
    maxMemory = 0;
  }
};

struct FootplanGoal {
  Vec2<float> goalPos;
};

using FootplanStateCost = float (*)(FootplanState&, FootplanGoal&);
using FootplanTransitionCost = float (*)(FootplanState&, FootplanState&, FootplanGoal&);

namespace FootplanCosts {
  float distanceToGoal(FootplanState& state, FootplanGoal& goal);
}

//  cheetah._bodyLength = 0.19 * 2;
//  cheetah._bodyWidth = 0.049 * 2;

class FootstepPlanner {
public:
  FootstepPlanner(bool verbose);
  void reset();
  void buildInputTrajectory(float duration, float dt, InputTrajectoryState x0, float omega);
  void planFixedEvenGait(std::vector<ContactState>& gait, float gait_period);
  std::vector<InputTrajectoryState>& getInitialTrajectory() {
    return _inputTrajectory;
  }

  void addCost(FootplanStateCost cost) {
    _stateCosts.push_back(cost);
  }

  void addCost(FootplanTransitionCost cost) {
    _transitionCosts.push_back(cost);
  }

  FootplanGoal& getGoal() {
    return _goal;
  }

  DefaultGaits defaults;
private:
  bool _verbose;
  FootplanStats _stats;
  FootplanGoal _goal;

  std::vector<FootplanStateCost> _stateCosts;
  std::vector<FootplanTransitionCost> _transitionCosts;
  std::vector<InputTrajectoryState> _inputTrajectory;
};


#endif //CHEETAH_SOFTWARE_GRAPHSEARCH_H
