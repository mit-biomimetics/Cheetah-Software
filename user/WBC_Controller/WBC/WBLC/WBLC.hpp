#ifndef WHOLE_BODY_LOCOMOTION_CONTROL_H
#define WHOLE_BODY_LOCOMOTION_CONTROL_H

#include <Utilities/Utilities_print.h>
#include <Goldfarb_Optimizer/QuadProg++.hh>
#include <WBC/ContactSpec.hpp>
#include <WBC/WBC.hpp>

template <typename T>
class WBLC_ExtraData {
 public:
  // Output
  DVec<T> opt_result_;
  DVec<T> qddot_;
  DVec<T> Fr_;

  // Input
  DVec<T> W_qddot_;
  DVec<T> W_rf_;
  DVec<T> W_xddot_;

  DVec<T> tau_min_;
  DVec<T> tau_max_;

  DVec<T> _des_jacc_cmd;

  WBLC_ExtraData() {}
  ~WBLC_ExtraData() {}
};

template <typename T>
class WBLC : public WBC<T> {
 public:
  WBLC(size_t num_qdot, const std::vector<ContactSpec<T>*>& contact_list);
  virtual ~WBLC() {}

  virtual void UpdateSetting(const DMat<T>& A, const DMat<T>& Ainv,
                             const DVec<T>& cori, const DVec<T>& grav,
                             void* extra_setting = NULL);

  virtual void MakeTorque(DVec<T>& cmd, void* extra_input = NULL);

 private:
  std::vector<ContactSpec<T>*> _contact_list;
  void _OptimizationPreparation(const DMat<T>& Aeq, const DVec<T>& beq,
                                const DMat<T>& Cieq, const DVec<T>& dieq);

  void _GetSolution(DVec<T>& cmd);
  void _OptimizationPreparation();

  size_t dim_opt_;
  size_t dim_eq_cstr_;     // equality constraints
  size_t dim_ieq_cstr_;    // inequality constraints
  size_t dim_first_task_;  // first task dimension
  WBLC_ExtraData<T>* data_;

  GolDIdnani::GVect<double> z;
  // Cost
  GolDIdnani::GMatr<double> G;
  GolDIdnani::GVect<double> g0;

  // Equality
  GolDIdnani::GMatr<double> CE;
  GolDIdnani::GVect<double> ce0;

  // Inequality
  GolDIdnani::GMatr<double> CI;
  GolDIdnani::GVect<double> ci0;

  size_t dim_rf_;
  size_t dim_Uf_;

  size_t dim_relaxed_task_;
  size_t dim_cam_;
  size_t dim_rf_cstr_;

  // BuildContactMtxVect builds the followings:
  void _BuildContactMtxVect();
  DMat<T> Uf_;
  DVec<T> Fr_ieq_;
  DVec<T> Fr_des_;
  DMat<T> Jc_;
  DVec<T> JcDotQdot_;

  // Setup the followings:
  void _Build_Equality_Constraint();
  DMat<T> Aeq_;
  DVec<T> beq_;

  void _Build_Inequality_Constraint();
  DMat<T> Cieq_;
  DVec<T> dieq_;

  DVec<T> qddot_;

  DMat<T> Sf_;  // floating base
  void _PrintDebug(T i) {
    (void)i;
    // printf("[WBLC] %f \n", i);
  }
};

#endif
