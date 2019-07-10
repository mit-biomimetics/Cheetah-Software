#include "FixedBodyContact.hpp"
#include <WBC_States/Cheetah_DynaCtrl_Definition.h>

template <typename T>
FixedBodyContact<T>::FixedBodyContact() : ContactSpec<T>(6) {
  Contact::Jc_ = DMat<T>::Zero(Contact::dim_contact_, cheetah::dim_config);

  for (size_t i(0); i < Contact::dim_contact_; ++i) Contact::Jc_(i, i) = 1.;
  Contact::Uf_ = DMat<T>::Zero(1, Contact::dim_contact_);
  Contact::ieq_vec_ = DVec<T>::Zero(1);
}

template <typename T>
FixedBodyContact<T>::~FixedBodyContact() {}

template <typename T>
bool FixedBodyContact<T>::_UpdateJc() {
  return true;
}

template <typename T>
bool FixedBodyContact<T>::_UpdateJcDotQdot() {
  Contact::JcDotQdot_ = DVec<T>::Zero(Contact::dim_contact_);
  return true;
}

template <typename T>
bool FixedBodyContact<T>::_UpdateUf() {
  return true;
}

template <typename T>
bool FixedBodyContact<T>::_UpdateInequalityVector() {
  return true;
}

template class FixedBodyContact<double>;
template class FixedBodyContact<float>;
