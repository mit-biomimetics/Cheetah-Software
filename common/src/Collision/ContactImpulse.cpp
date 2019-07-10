#include "Collision/ContactImpulse.h"
#include "Utilities/Utilities_print.h"

template <typename T>
void ContactImpulse<T>::UpdateQdot(FBModelState<T>& state) {
  CC::_nContact = CC::_CheckContact();
  if (CC::_nContact > 0) {
    DVec<T> qdot(_nDof);

    for (size_t i(0); i < 6; ++i) {
      qdot[i] = state.bodyVelocity[i];
    }
    for (size_t i(0); i < _nDof - 6; ++i) {
      qdot[i + 6] = state.qd[i];
    }

    _UpdateVelocity(qdot);

    for (size_t i(0); i < 6; ++i) {
      state.bodyVelocity[i] = qdot[i];
    }
    for (size_t i(0); i < _nDof - 6; ++i) {
      state.qd[i] = qdot[i + 6];
    }

    // Update contact force w.r.t. global frame
    // This is not neccessary for dynamics, but
    // the global contact forces are used for
    // graphics
    for (size_t i(0); i < CC::_nContact; ++i) {
      // Divide by dt to convert impulses to forces
      CC::_cp_force_list[CC::_idx_list[i]] =
          CC::_cp_frame_list[i] * CC::_cp_local_force_list[i] / _dt;
    }
  }
}

template <typename T>
void ContactImpulse<T>::_UpdateVelocity(DVec<T>& qdot) {
  T* lambda_list_x = new T[CC::_nContact];
  T* lambda_list_y = new T[CC::_nContact];
  T* lambda_list_z = new T[CC::_nContact];

  T* des_vel_list_z = new T[CC::_nContact];
  T* des_vel_list_tan = new T[CC::_nContact];

  T* min_list_z = new T[CC::_nContact];
  T* max_list_z = new T[CC::_nContact];

  T* min_list_tan = new T[CC::_nContact];
  T* max_list_tan = new T[CC::_nContact];

  vectorAligned<DVec<T> > AinvB_list_x(CC::_nContact);
  vectorAligned<DVec<T> > AinvB_list_y(CC::_nContact);
  vectorAligned<DVec<T> > AinvB_list_z(CC::_nContact);

  vectorAligned<D3Mat<T> > Jc_list(CC::_nContact);
  // For check
  vectorAligned<Vec3<T> > cp_vel_list(CC::_nContact);

  DVec<T> AinvB;
  Vec3<T> cp_local_vel;
  Vec3<T> direction;
  CC::_model->contactJacobians();

  // Prepare Matrix and Vector
  for (size_t i(0); i < CC::_nContact; ++i) {
    // Lambda & Ainv*J'
    lambda_list_x[i] = CC::_model->applyTestForce(
        CC::_idx_list[i], CC::_cp_frame_list[i].template block<3, 1>(0, 0),
        AinvB_list_x[i]);
    lambda_list_x[i] = 1. / lambda_list_x[i];

    lambda_list_y[i] = CC::_model->applyTestForce(
        CC::_idx_list[i], CC::_cp_frame_list[i].template block<3, 1>(0, 1),
        AinvB_list_y[i]);
    lambda_list_y[i] = 1. / lambda_list_y[i];

    lambda_list_z[i] = CC::_model->applyTestForce(
        CC::_idx_list[i], CC::_cp_frame_list[i].template block<3, 1>(0, 2),
        AinvB_list_z[i]);
    lambda_list_z[i] = 1. / lambda_list_z[i];

    // Local Velocity
    cp_local_vel =
        CC::_cp_frame_list[i].transpose() * CC::_model->_vGC[CC::_idx_list[i]];

    // Desired Velocity
    des_vel_list_tan[i] = 0.;
    if (cp_local_vel[2] < 0.) {
      des_vel_list_z[i] =
          -CC::_cp_resti_list[i] * cp_local_vel[2] -
          _penetration_recover_ratio * CC::_cp_penetration_list[i];
    } else {
      des_vel_list_z[i] =
          std::max(cp_local_vel[2],
                   -_penetration_recover_ratio * CC::_cp_penetration_list[i]);
    }

    // Contact Jacobian
    Jc_list[i] =
        CC::_cp_frame_list[i].transpose() * CC::_model->_Jc[CC::_idx_list[i]];

    min_list_z[i] = 0.;
    max_list_z[i] = 1.e5;
  }

  // Update Velocity & Find Impulse Force
  for (size_t iter(0); iter < _iter_lim; ++iter) {
    _iter_sum = 0;
    // Normal (Z) *********************************************************
    _UpdateQdotOneDirection(2, Jc_list, lambda_list_z, AinvB_list_z,
                            des_vel_list_z, min_list_z, max_list_z, qdot);

    // X ******************************************************************
    for (size_t i(0); i < CC::_nContact; ++i) {
      max_list_tan[i] = CC::_cp_mu_list[i] * CC::_cp_local_force_list[i][2];
      min_list_tan[i] = -max_list_tan[i];
    }
    _UpdateQdotOneDirection(0, Jc_list, lambda_list_x, AinvB_list_x,
                            des_vel_list_tan, min_list_tan, max_list_tan, qdot);

    // Y ******************************************************************
    _UpdateQdotOneDirection(1, Jc_list, lambda_list_y, AinvB_list_y,
                            des_vel_list_tan, min_list_tan, max_list_tan, qdot);

    if (_iter_sum < 5) {
      // printf("converged: %lu \n", _iter_sum);
      break;
    }
    _iter_sum = 0;
  }

  delete[] lambda_list_x;
  delete[] lambda_list_y;
  delete[] lambda_list_z;

  delete[] des_vel_list_z;
  delete[] des_vel_list_tan;

  delete[] min_list_z;
  delete[] max_list_z;

  delete[] min_list_tan;
  delete[] max_list_tan;
}

template <typename T>
void ContactImpulse<T>::_UpdateQdotOneDirection(
    size_t idx, const vectorAligned<D3Mat<T> >& Jc_list, const T* lambda_list,
    const vectorAligned<DVec<T> > AinvB_list, const T* des_vel_list,
    const T* min_list, const T* max_list, DVec<T>& qdot) {
  T dforce, pre_force, cp_vel, dforce_sum;
  for (size_t iter(0); iter < _iter_lim; ++iter) {
    dforce_sum = 0;
    for (size_t i(0); i < CC::_nContact; ++i) {
      cp_vel = (Jc_list[i].block(idx, 0, 1, _nDof) * qdot)(0, 0);

      dforce = (des_vel_list[i] - cp_vel) * lambda_list[i];
      pre_force = CC::_cp_local_force_list[i][idx];
      CC::_cp_local_force_list[i][idx] =
          std::max(min_list[i], std::min(pre_force + dforce, max_list[i]));

      dforce = CC::_cp_local_force_list[i][idx] - pre_force;

      qdot += (AinvB_list[i] * dforce);
      dforce_sum += fabs(dforce);
    }
    if (dforce_sum < _tol) {
      _iter_sum += iter;
      // printf("dforce sum is small enough: %f\n", dforce_sum);
      break;
    }
  }
}

template class ContactImpulse<double>;
template class ContactImpulse<float>;
