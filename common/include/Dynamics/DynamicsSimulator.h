/*! @file DynamicsSimulator.h
 *  @brief Rigid Body Dynamics Simulator with Collisions
 *
 *  Combines ABA, Collisions, integrator, and any other external forces to run a
 * simulation. Doesn't do any graphics.
 */

#ifndef PROJECT_DYNAMICSSIMULATOR_H
#define PROJECT_DYNAMICSSIMULATOR_H

#include "cppTypes.h"
#include "Collision/CollisionBox.h"
#include "Collision/CollisionMesh.h"
#include "Collision/CollisionPlane.h"
#include "Collision/ContactConstraint.h"
#include "FloatingBaseModel.h"
#include "Math/orientation_tools.h"
#include "spatial.h"

using namespace ori;
using namespace spatial;

/*!
 * Parameters used to move the robot to a "home" position
 */
template <typename T>
struct RobotHomingInfo {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vec3<T> position;
  Vec3<T> rpy;  // body coordinates
  double kp_lin;
  double kd_lin;
  double kp_ang;
  double kd_ang;
  bool active_flag;
};

/*!
 * Class (containing state) for dynamics simulation of a floating-base system
 */
template <typename T>
class DynamicsSimulator {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  DynamicsSimulator(
      FloatingBaseModel<T>& model,
      bool useSpringDamper = false);  //! Initialize simulator with given model
  void step(T dt, const DVec<T>& tau, T kp,
            T kd);  //! Simulate forward one step

  //! Find _dstate with the articulated body algorithm
  void runABA(const DVec<T>& tau) { _model.runABA(tau, _dstate); }

  //! Do forward kinematics for feet
  void forwardKinematics() { _model.forwardKinematics(); }

  //! Integrate to find new _state
  void integrate(T dt);

  /*!
   * Set the state of the robot being simulated
   */
  void setState(const FBModelState<T>& state) {
    _state = state;
    _model.setState(state);  // force recalculate dynamics
  }

  /*!
   * Set the homing parameters of the simulator
   * @param homing : homing parameters
   */
  void setHoming(const RobotHomingInfo<T>& homing) {
    _homing = homing;
  }

  /*!
   * Get the state of the robot
   * @return
   */
  const FBModelState<T>& getState() const { return _state; }

  /*!
   * Get the most recently calculated state derivative (updated by runABA)
   * @return
   */
  const FBModelStateDerivative<T>& getDState() const { return _dstate; }

  /*!
   * Add external forces. These are added on top of the ground reaction forces
   * The i-th force is the spatial force applied to body i.
   * This is cleared after each step of the simulator.
   */
  void setAllExternalForces(const vectorAligned<SVec<T>>& forces) {
    _model._externalForces = forces;
  }

  /*!
   * Add a collision plane (xy plane)
   * @param mu : friction coefficient
   * @param rest : restitution coefficient
   * @param height : height of plane
   */
  void addCollisionPlane(T mu, T rest, T height) {
    _contact_constr->AddCollision(new CollisionPlane<T>(mu, rest, height));
  }

  /*!
   * Add a box
   * @param mu : friction coefficient
   * @param rest : restitution coefficient
   * @param depth : depth of box
   * @param width : width of box
   * @param height : height of box
   * @param pos : position of center of box
   * @param ori : orientation of box (rotation matrix)
   */
  void addCollisionBox(T mu, T rest, T depth, T width, T height,
                       const Vec3<T>& pos, const Mat3<T>& ori) {
    _contact_constr->AddCollision(
        new CollisionBox<T>(mu, rest, depth, width, height, pos, ori));
  }

  /*!
   * Add a height map mesh
   * @param mu : friction coefficient
   * @param rest : restitution coefficient
   * @param grid_size : size of each grid (must be square)
   * @param left_corner_loc : location of 0,0 coordinate
   * @param height_map : height map data
   */
  void addCollisionMesh(T mu, T rest, T grid_size,
                        const Vec3<T>& left_corner_loc,
                        const DMat<T>& height_map) {
    _contact_constr->AddCollision(
        new CollisionMesh<T>(mu, rest, grid_size, left_corner_loc, height_map));
  }

  /*!
   * Get the number of bodies (not including rotors)
   * @return number of bodies
   */
  size_t getNumBodies() { return _model._nDof; }

  /*!
   * Get the number of maximum possible ground contact points
   * @return The number of maximum possible ground contact points
   */
  const size_t& getTotalNumGC() { return _model._nGroundContact; }

  /*!
   * Get the force at a given contact point
   * @param idx : The index of the contact point
   * @return The force at the contact point (global frame)
   */
  const Vec3<T>& getContactForce(size_t idx) {
    return _contact_constr->getGCForce(idx);
  }

  /*!
   * Get the floating base model which is being simulated
   * @return floating base model
   */
  const FloatingBaseModel<T>& getModel() { return _model; }

 private:
  void updateCollisions(T dt, T kp, T kd);  //! Update ground collision list
  FBModelState<T> _state;
  FBModelStateDerivative<T> _dstate;
  FloatingBaseModel<T>& _model;
  vector<CollisionPlane<T>> _collisionPlanes;
  ContactConstraint<T>* _contact_constr;
  SVec<T> _lastBodyVelocity;
  bool _useSpringDamper;

  RobotHomingInfo<T> _homing;
};

#endif  // PROJECT_DYNAMICSSIMULATOR_H
