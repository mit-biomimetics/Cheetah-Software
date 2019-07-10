#include "SingleContact.hpp"
#include <Utilities/Utilities_print.h>
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

// [ Fx, Fy, Fz ]
template <typename T>
SingleContact<T>::SingleContact(const FloatingBaseModel<T>* robot, int pt)
    : ContactSpec<T>(3), _max_Fz(1500.), _contact_pt(pt), _dim_U(6) {
  Contact::idx_Fz_ = 2;
  robot_sys_ = robot;
  Contact::Jc_ = DMat<T>(Contact::dim_contact_, cheetah::dim_config);
  Contact::JcDotQdot_ = DVec<T>::Zero(Contact::dim_contact_);
  Contact::Uf_ = DMat<T>::Zero(_dim_U, Contact::dim_contact_);

  T mu(0.3);

  Contact::Uf_(0, 2) = 1.;

  Contact::Uf_(1, 0) = 1.;
  Contact::Uf_(1, 2) = mu;
  Contact::Uf_(2, 0) = -1.;
  Contact::Uf_(2, 2) = mu;

  Contact::Uf_(3, 1) = 1.;
  Contact::Uf_(3, 2) = mu;
  Contact::Uf_(4, 1) = -1.;
  Contact::Uf_(4, 2) = mu;

  // Upper bound of normal force
  Contact::Uf_(5, 2) = -1.;
}

template <typename T>
SingleContact<T>::~SingleContact() {}

template <typename T>
bool SingleContact<T>::_UpdateJc() {
  Contact::Jc_ = robot_sys_->_Jc[_contact_pt];

  // Quat<T> quat = robot_sys_->_state.bodyOrientation;
  // Mat3<T> Rot = ori::quaternionToRotationMatrix(quat);
  // Contact::Jc_.block(0,3, 3,3) = Rot*Contact::Jc_.block(0,3,3,3);

  // Contact::Jc_.block(0,0, 3,3) = Rot.transpose()*Contact::Jc_.block(0,0,3,3);
  // pretty_print(Rot, std::cout, "body ori");
  // pretty_print(Contact::Jc_, std::cout, "Jc");
  return true;
}

template <typename T>
bool SingleContact<T>::_UpdateJcDotQdot() {
  Contact::JcDotQdot_ = robot_sys_->_Jcdqd[_contact_pt];
  // pretty_print(Contact::JcDotQdot_, std::cout, "JcDotQdot");
  return true;
}

template <typename T>
bool SingleContact<T>::_UpdateUf() {
  return true;
}

template <typename T>
bool SingleContact<T>::_UpdateInequalityVector() {
  Contact::ieq_vec_ = DVec<T>::Zero(_dim_U);
  Contact::ieq_vec_[5] = -_max_Fz;
  return true;
}

template class SingleContact<double>;
template class SingleContact<float>;
