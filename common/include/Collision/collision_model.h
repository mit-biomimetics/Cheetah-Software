/*! @file collision_model.h
 *  @brief Nonlinear friction and collision model used in the simulator for ground contact
 *
 *  The collision model is based on the gcontact.m file in spatial_v2 and modified by Pat for codegen
 *  This uses the model in "Modelling the Contact Between a Rolling Sphere and a Compliant Ground Plane" by
 *  Azad and Featherstone.  It assumes a point contact.
 *
 */

#ifndef PROJECT_COLLISION_MODEL_H
#define PROJECT_COLLISION_MODEL_H

#include "cppTypes.h"
#include "Dynamics/spatial.h"
#include <cmath>

using namespace spatial;


/*!
 * Run the ground contact model for a single collision plane on a list of ground contact points
 * The ground is allowed to deform in the tangential direction, but not the normal direction.
 * The ground also "remembers" its deformation between separate contact events. (however it does spring back pretty quickly)
 * @param _pGC List of ground contact point locations.
 * @param _vGC List of ground contact point velocities
 * @param _deflections List of deflection state variables.  This is updated by this model
 * @param K Ground stiffness
 * @param D Ground damping
 * @param mu Ground friction
 * @param dt Timestep (used for deflection)
 * @param X The coordinate transformation to the ground plane 
 * (the ground plane is the xy-plane in its coordinates)
 */
template <typename T>
void groundContactModelWithOffset(vectorAligned<Vec3<T>>& _pGC, vectorAligned<Vec3<T>>& _vGC,
                                  vectorAligned<Vec2<T>>& _deflections, 
                                  vectorAligned<Vec3<T>>& _forces, T K, T D, T mu, T dt, Mat6<T>& X) {

  Mat3<T> R = rotationFromSXform(X); // rotation to ground plane
  size_t nPt = _pGC.size();          // number of ground contact points
  for(size_t i = 0; i < nPt; i++) {
    Vec3<T> p = sXFormPoint(X, _pGC[i]);  // position in plane coordinates
    Vec3<T> v = R * _vGC[i];              // velocity in plane coordinates
    Vec2<T> deflection = _deflections[i]; // the deflection of the ground from previously
    T z = p[2];                           // the penetration into the ground
    T zd = v[2];                          // the penetration velocity
    T zr = std::sqrt(std::max(T(0), -z)); // sqrt penetration into the ground, or zero if we aren't penetrating
    T normalForce = zr * (-K*z - D*zd);   // normal force is spring-damper * sqrt(penetration)
    bool inContact = normalForce > 0;     // contact flag

    // set output force to zero for now.
    _forces[i][0] = 0;
    _forces[i][1] = 0;
    _forces[i][2] = 0;

    // first assume there's no contact, so the ground "springs back"
    Vec2<T> deflectionRate = (-K/D) * deflection.template topLeftCorner<2,1>();

    if(inContact) {
      _forces[i][2] = normalForce;   // set the normal force. This is in the plane's coordinates for now

      // first, assume sticking
      // this means the tangential deformation happens at the speed of the foot.
      deflectionRate = v.template topLeftCorner<2,1>();
      Vec2<T> tangentialSpringForce = K * zr * deflection; // tangential force due to "spring"
      Vec2<T> tangentialForce = -tangentialSpringForce - D * zr * deflectionRate; // add damping to get total tangential

      // check for slipping:
      T slipForce = mu * normalForce;        // maximum force magnitude without slipping
      T tangentialForceMagnitude = tangentialForce.norm(); // actual force magnitude if we assume sticking
      T r = tangentialForceMagnitude / slipForce;          // ratio of force/max_force

      if(r > 1) {
        // we are slipping.
        tangentialForce = tangentialForce / r;    // adjust tangential force to avoid slipping
        deflectionRate = - (tangentialForce + tangentialSpringForce) / (D * zr);
      }
      // set forces
      _forces[i][0] = tangentialForce[0];
      _forces[i][1] = tangentialForce[1];
    }
    // integrate ground deflection
    _deflections[i] += dt * deflectionRate;
    // move back into robot frame
    _forces[i] = R.transpose() * _forces[i];
  }
}


/*!
 * Run the ground contact model on a series of points.
 * This assumes that the ground is the xy-plane.
 * The ground is allowed the deform in the tangential direction, but not the normal direction.
 * The ground also "remembers" its deformation between separate contact events. (however it does spring back pretty quickly)
 * @param _pGC List of ground contact point locations.
 * @param _vGC List of ground contact point velocities
 * @param _deflections List of deflection state variables.  This is updated by this model
 * @param K Ground stiffness
 * @param D Ground damping
 * @param mu Ground friction
 * @param dt Timestep (used for deflection)
 */
template <typename T>
void groundContactModel(vectorAligned<Vec3<T>>& _pGC, vectorAligned<Vec3<T>>& _vGC,
        vectorAligned<Vec2<T>>& _deflections, vectorAligned<Vec3<T>>& _forces, T K, T D, T mu, T dt) {
  Mat6<T> eye = Mat6<T>::Identity();
  groundContactModelWithOffset(_pGC, _vGC, _deflections, _forces, K, D, mu, dt, eye);
}





#endif //PROJECT_COLLISION_MODEL_H
