#ifndef CONTACT_SPEC
#define CONTACT_SPEC

#include <cppTypes.h>

#define Contact ContactSpec<T>

template <typename T>
class ContactSpec {
 public:
  ContactSpec(size_t dim) : dim_contact_(dim), b_set_contact_(false) {
    idx_Fz_ = dim - 1;  // because normally (tau_x,y,z , linear_x,y,z)
    Fr_des_ = DVec<T>::Zero(dim);
  }
  virtual ~ContactSpec() {}

  size_t getDim() const { return dim_contact_; }
  size_t getDimRFConstraint() const { return Uf_.rows(); }
  size_t getFzIndex() const { return idx_Fz_; }

  void getContactJacobian(DMat<T>& Jc) { Jc = Jc_; }
  void getJcDotQdot(DVec<T>& JcDotQdot) { JcDotQdot = JcDotQdot_; }
  void UnsetContact() { b_set_contact_ = false; }

  void getRFConstraintMtx(DMat<T>& Uf) { Uf = Uf_; }
  void getRFConstraintVec(DVec<T>& ieq_vec) { ieq_vec = ieq_vec_; }
  const DVec<T>& getRFDesired() { return Fr_des_; }
  void setRFDesired(const DVec<T>& Fr_des) { Fr_des_ = Fr_des; }

  bool UpdateContactSpec() {
    _UpdateJc();
    _UpdateJcDotQdot();
    _UpdateUf();
    _UpdateInequalityVector();
    b_set_contact_ = true;
    return true;
  }

 protected:
  virtual bool _UpdateJc() = 0;
  virtual bool _UpdateJcDotQdot() = 0;
  virtual bool _UpdateUf() = 0;
  virtual bool _UpdateInequalityVector() = 0;

  int idx_Fz_;
  DMat<T> Uf_;
  DVec<T> ieq_vec_;
  DVec<T> Fr_des_;

  DMat<T> Jc_;
  DVec<T> JcDotQdot_;
  size_t dim_contact_;
  bool b_set_contact_;
};
#endif
