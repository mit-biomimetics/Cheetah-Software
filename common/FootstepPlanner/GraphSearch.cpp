#include "GraphSearch.h"
#include "Math/orientation_tools.h"

float FootplanCosts::distanceToGoal(FootplanState &state, FootplanGoal &goal) {
  Vec2<float> dp = state.pBase - goal.goalPos;
  return dp.norm();
}

FootstepPlanner::FootstepPlanner(bool verbose) : _verbose(verbose) {
  _stats.reset();
  defaults.trotting = {{true, false, false, true},
                       {false, true, true, false}};
  defaults.standing = {{true, true, true, true}};
}

void FootstepPlanner::reset() {
  _stats.reset();
  _stateCosts.clear();
  _transitionCosts.clear();
}

void FootstepPlanner::buildInputTrajectory(float duration, float dt, InputTrajectoryState x0, float omega) {
  if(_verbose) {
    printf("Input trajectory with %d steps\n", (int)(duration / dt));
  }
  _inputTrajectory.clear();
  _inputTrajectory.reserve(duration / dt);

  Vec3<float> velocity(x0.v[0], x0.v[1], 0.f);
  Vec3<float> position(x0.p[0], x0.p[1], 0.f);
  float theta = x0.theta;
  float t = 0;
  for(uint32_t i = 0; i < (duration / dt); i++) {
    Vec3<float> vRot = ori::coordinateRotation(ori::CoordinateAxis::Z, theta).transpose() * velocity;

    _inputTrajectory.push_back({{position[0], position[1]}, {vRot[0], vRot[1]}, theta});

    position += vRot * dt;
    t += dt;
    theta += omega * dt;
  }
}

void FootstepPlanner::planFixedEvenGait(std::vector<ContactState> &gait, float gait_period) {
  (void)gait;
  (void)gait_period;
}