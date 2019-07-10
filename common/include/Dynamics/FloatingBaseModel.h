/*! @file FloatingBaseModel.h
 *  @brief Implementation of Rigid Body Floating Base model data structure
 *
 * This class stores the kinematic tree described in "Rigid Body Dynamics
 * Algorithms" by Featherstone (download from
 * https://www.springer.com/us/book/9780387743141 on MIT internet)
 *
 * The tree includes an additional "rotor" body for each body.  This rotor is
 * fixed to the parent body and has a gearing constraint.  This is efficiently
 * included using a technique similar to what is described in Chapter 12 of
 * "Robot and Multibody Dynamics" by Jain.  Note that this implementation is
 * highly specific to the case of a single rotating rotor per rigid body. Rotors
 * have the same joint type as their body, but with an additional gear ratio
 * multiplier applied to the motion subspace. The rotors associated with the
 * floating base don't do anything.
 */

#ifndef LIBBIOMIMETICS_FLOATINGBASEMODEL_H
#define LIBBIOMIMETICS_FLOATINGBASEMODEL_H

#include "Math/orientation_tools.h"
#include "SpatialInertia.h"
#include "spatial.h"

#include <eigen3/Eigen/StdVector>

#include <string>
#include <vector>

using std::vector;
using namespace ori;
using namespace spatial;

/*!
 * The state of a floating base model (base and joints)
 */
template <typename T>
struct FBModelState {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Quat<T> bodyOrientation;
  Vec3<T> bodyPosition;
  SVec<T> bodyVelocity;  // body coordinates
  DVec<T> q;
  DVec<T> qd;

  void print() const {
    printf("position: %.3f %.3f %.3f\n", bodyPosition[0], bodyPosition[1],
           bodyPosition[2]);
  }
};

/*!
 * The result of running the articulated body algorithm on a rigid-body floating
 * base model
 */
template <typename T>
struct FBModelStateDerivative {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vec3<T> dBodyPosition;
  SVec<T> dBodyVelocity;
  DVec<T> qdd;
};

/*!
 * Class to represent a floating base rigid body model with rotors and ground
 * contacts. No concept of state.
 */
template <typename T>
class FloatingBaseModel {
 public:
  /*!
   * Initialize a floating base model with default gravity
   */
  FloatingBaseModel() : _gravity(0, 0, -9.81) {}
  ~FloatingBaseModel() {}

  /*!
   * Add floating base.  Must be the first body added, and there can only be one
   */
  void addBase(const SpatialInertia<T>& inertia);

  /*!
   * Add floating base.  Must be the first body added, and there can only be one
   */
  void addBase(T mass, const Vec3<T>& com, const Mat3<T>& I);

  /*!
   * Add a point for collisions
   * @param bodyID   : body that the point belongs to (body 5 for floating base)
   * @param location : location of point in body coordinates
   * @param isFoot   : if the point is a foot or not.  Only feet have their
   * Jacobian calculated on the robot
   * @return collisionPointID of the new point
   */
  int addGroundContactPoint(int bodyID, const Vec3<T>& location,
                            bool isFoot = false);

  /*!
   * Add bounding box collision points around a body.
   * @param bodyId : Body to add
   * @param dims   : Dimension of points
   */
  void addGroundContactBoxPoints(int bodyId, const Vec3<T>& dims);

  /*!
   * Add a body to the tree
   * @param inertia        : Inertia of body (body coords)
   * @param rotorInertia   : Inertia of rotor (rotor coords)
   * @param gearRatio      : Gear ratio.  >1 for a gear reduction
   * @param parent         : Body ID of the link that the body is connected to
   * @param jointType      : Type of joint
   * @param jointAxis      : Axis of joint
   * @param Xtree          : Location of joint
   * @param Xrot           : Location of rotor
   * @return               : bodyID
   */
  int addBody(const SpatialInertia<T>& inertia,
              const SpatialInertia<T>& rotorInertia, T gearRatio, int parent,
              JointType jointType, CoordinateAxis jointAxis,
              const Mat6<T>& Xtree, const Mat6<T>& Xrot);

  /*!
   * Add a body to the tree
   * @param inertia        : Inertia of body (body coords)
   * @param rotorInertia   : Inertia of rotor (rotor coords)
   * @param gearRatio      : Gear ratio.  >1 for a gear reduction
   * @param parent         : Body ID of the link that the body is connected to
   * @param jointType      : Type of joint
   * @param jointAxis      : Axis of joint
   * @param Xtree          : Location of joint
   * @param Xrot           : Location of rotor
   * @return               : bodyID
   */
  int addBody(const MassProperties<T>& inertia,
              const MassProperties<T>& rotorInertia, T gearRatio, int parent,
              JointType jointType, CoordinateAxis jointAxis,
              const Mat6<T>& Xtree, const Mat6<T>& Xrot);

  /*!
   * Very simple to check to make sure the dimensions are correct
   */
  void check();

  /*!
   * Total mass of all rotors
   */
  T totalRotorMass();

  /*!
   * Total mass of all bodies which are not rotors
   */
  T totalNonRotorMass();

  const std::vector<int>& getParentVector() { return _parents; }

  const std::vector<SpatialInertia<T>,
                    Eigen::aligned_allocator<SpatialInertia<T>>>&
  getBodyInertiaVector() {
    return _Ibody;
  }

  const std::vector<SpatialInertia<T>,
                    Eigen::aligned_allocator<SpatialInertia<T>>>&
  getRotorInertiaVector() {
    return _Irot;
  }

  /*!
   * Set the gravity
   */
  void setGravity(Vec3<T>& g) { _gravity = g; }

  void setContactComputeFlag(size_t gc_index, bool flag) {
    _compute_contact_info[gc_index] = flag;
  }

  DMat<T> invContactInertia(const int gc_index,
                            const D6Mat<T>& force_directions);
  T invContactInertia(const int gc_index, const Vec3<T>& force_ics_at_contact);

  T applyTestForce(const int gc_index, const Vec3<T>& force_ics_at_contact,
                   FBModelStateDerivative<T>& dstate_out);

  T applyTestForce(const int gc_index, const Vec3<T>& force_ics_at_contact,
                   DVec<T>& dstate_out);

  void addDynamicsVars(int count);

  void resizeSystemMatricies();

  void setState(const FBModelState<T>& state) {
    _state = state;

    _biasAccelerationsUpToDate = false;
    _compositeInertiasUpToDate = false;

    resetCalculationFlags();
  }
  void resetCalculationFlags() {
    _articulatedBodiesUpToDate = false;
    _kinematicsUpToDate = false;
    _forcePropagatorsUpToDate = false;
    _qddEffectsUpToDate = false;
    _accelerationsUpToDate = false;
  }

  void setDState(const FBModelStateDerivative<T>& dState) {
    _dState = dState;
    _accelerationsUpToDate = false;
  }

  Vec3<T> getPosition(const int link_idx, const Vec3<T> & local_pos);
  Vec3<T> getPosition(const int link_idx);


  Mat3<T> getOrientation(const int link_idx);
  Vec3<T> getLinearVelocity(const int link_idx, const Vec3<T>& point);
  Vec3<T> getLinearVelocity(const int link_idx);

  Vec3<T> getLinearAcceleration(const int link_idx, const Vec3<T>& point);
  Vec3<T> getLinearAcceleration(const int link_idx);

  Vec3<T> getAngularVelocity(const int link_idx);
  Vec3<T> getAngularAcceleration(const int link_idx);

  void forwardKinematics();
  void biasAccelerations();
  void compositeInertias();
  void forwardAccelerationKinematics();
  void contactJacobians();

  DVec<T> generalizedGravityForce();
  DVec<T> generalizedCoriolisForce();
  DMat<T> massMatrix();
  DVec<T> inverseDynamics(const FBModelStateDerivative<T>& dState);
  void runABA(const DVec<T>& tau, FBModelStateDerivative<T>& dstate);

  size_t _nDof = 0;
  Vec3<T> _gravity;
  vector<int> _parents;
  vector<T> _gearRatios;
  vector<T> _d, _u;

  vector<JointType> _jointTypes;
  vector<CoordinateAxis> _jointAxes;
  vector<Mat6<T>, Eigen::aligned_allocator<Mat6<T>>> _Xtree, _Xrot;
  vector<SpatialInertia<T>, Eigen::aligned_allocator<SpatialInertia<T>>> _Ibody,
      _Irot;
  vector<std::string> _bodyNames;

  size_t _nGroundContact = 0;
  vector<size_t> _gcParent;
  vector<Vec3<T>> _gcLocation;
  vector<uint64_t> _footIndicesGC;

  vector<Vec3<T>> _pGC;
  vector<Vec3<T>> _vGC;

  vector<bool> _compute_contact_info;

  const DMat<T>& getMassMatrix() const { return _H; }
  const DVec<T>& getGravityForce() const { return _G; }
  const DVec<T>& getCoriolisForce() const { return _Cqd; }

  // void getPositionVelocity(
  // const int & link_idx, const Vec3<T> & local_pos,
  // Vec3<T> & link_pos, Vec3<T> & link_vel) const ;

  /// BEGIN ALGORITHM SUPPORT VARIABLES
  FBModelState<T> _state;
  FBModelStateDerivative<T> _dState;

  vectorAligned<SVec<T>> _v, _vrot, _a, _arot, _avp, _avprot, _c, _crot, _S,
      _Srot, _fvp, _fvprot, _ag, _agrot, _f, _frot;

  vectorAligned<SVec<T>> _U, _Urot, _Utot, _pA, _pArot;
  vectorAligned<SVec<T>> _externalForces;

  vectorAligned<SpatialInertia<T>> _IC;
  vectorAligned<Mat6<T>> _Xup, _Xa, _Xuprot, _IA, _ChiUp;

  DMat<T> _H, _C;
  DVec<T> _Cqd, _G;

  vectorAligned<D6Mat<T>> _J;
  vectorAligned<SVec<T>> _Jdqd;

  vectorAligned<D3Mat<T>> _Jc;
  vectorAligned<Vec3<T>> _Jcdqd;

  bool _kinematicsUpToDate = false;
  bool _biasAccelerationsUpToDate = false;
  bool _accelerationsUpToDate = false;

  bool _compositeInertiasUpToDate = false;

  void updateArticulatedBodies();
  void updateForcePropagators();
  void udpateQddEffects();

  void resetExternalForces() {
    for (size_t i = 0; i < _nDof; i++) {
      _externalForces[i] = SVec<T>::Zero();
    }
  }

  bool _articulatedBodiesUpToDate = false;
  bool _forcePropagatorsUpToDate = false;
  bool _qddEffectsUpToDate = false;

  DMat<T> _qdd_from_base_accel;
  DMat<T> _qdd_from_subqdd;
  Eigen::ColPivHouseholderQR<Mat6<T>> _invIA5;
};

#endif  // LIBBIOMIMETICS_FLOATINGBASEMODEL_H
