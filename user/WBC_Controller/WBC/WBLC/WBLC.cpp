#include "WBLC.hpp"
#include <Utilities/Timer.h>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>

#define _WBLC_TORQUE_LIMIT_ false

template <typename T>
WBLC<T>::WBLC(size_t num_qdot, const std::vector<ContactSpec<T>*>& contact_list)
    : WBC<T>(num_qdot) {
  dim_rf_ = 0;
  dim_Uf_ = 0;
  for (size_t i(0); i < (contact_list).size(); ++i) {
    _contact_list.push_back(contact_list[i]);
    dim_rf_ += (_contact_list)[i]->getDim();
    dim_Uf_ += (_contact_list)[i]->getDimRFConstraint();
  }
  // Dimension Setting
  dim_eq_cstr_ = 6 + dim_rf_;
#if (_WBLC_TORQUE_LIMIT_)
  dim_opt_ = WB::num_qdot_ + 2 * dim_rf_;  // (delta_qddot, Fr, xddot_c)
  dim_ieq_cstr_ =
      2 * WB::num_act_joint_ + dim_Uf_;  // torque limit, friction cone
#else
  dim_opt_ = 6 + 2 * dim_rf_;  // (delta_virtual, Fr, xddot_c)
  dim_ieq_cstr_ = dim_Uf_;     // friction cone
#endif
  // dim_ieq_cstr_ = dim_Uf_; // torque limit, friction cone
  // printf("dim: opt, eq, ieq = %lu, %lu, %lu\n", dim_opt_, dim_eq_cstr_,
  // dim_ieq_cstr_);

  qddot_ = DVec<T>::Zero(WB::num_qdot_);

  G.resize(dim_opt_, dim_opt_);
  g0.resize(dim_opt_);
  CE.resize(dim_opt_, dim_eq_cstr_);
  ce0.resize(dim_eq_cstr_);
  CI.resize(dim_opt_, dim_ieq_cstr_);
  ci0.resize(dim_ieq_cstr_);

  Aeq_ = DMat<T>::Zero(dim_eq_cstr_, dim_opt_);
  beq_ = DVec<T>::Zero(dim_eq_cstr_);
  Cieq_ = DMat<T>::Zero(dim_ieq_cstr_, dim_opt_);
  dieq_ = DVec<T>::Zero(dim_ieq_cstr_);

  Jc_ = DMat<T>::Zero(dim_rf_, WB::num_qdot_);
  JcDotQdot_ = DVec<T>::Zero(dim_rf_);
  Fr_des_ = DVec<T>::Zero(dim_rf_);
  Uf_ = DMat<T>::Zero(dim_Uf_, dim_rf_);
  Fr_ieq_ = DVec<T>::Zero(dim_Uf_);

  for (size_t i(0); i < dim_opt_; ++i) {
    for (size_t j(0); j < dim_opt_; ++j) {
      G[i][j] = 0.;
    }
    g0[i] = 0.;
  }
}

template <typename T>
void WBLC<T>::UpdateSetting(const DMat<T>& A, const DMat<T>& Ainv,
                            const DVec<T>& cori, const DVec<T>& grav,
                            void* extra_setting) {
  WB::A_ = A;
  WB::Ainv_ = Ainv;
  WB::cori_ = cori;
  WB::grav_ = grav;
  WB::b_updatesetting_ = true;

  (void)extra_setting;
  // pretty_print(grav_, std::cout, "grav");
  // pretty_print(cori_, std::cout, "cori");
  // pretty_print(WB::A_, std::cout, "A");
}

template <typename T>
void WBLC<T>::MakeTorque(DVec<T>& cmd, void* extra_input) {
  if (!WB::b_updatesetting_) {
    printf("[Wanning] WBLC setting is not done\n");
  }
  if (extra_input) data_ = static_cast<WBLC_ExtraData<T>*>(extra_input);

  qddot_.setZero();
  for (size_t i(0); i < WB::num_act_joint_; ++i) {
    //(void)des_jacc_cmd[i];
    qddot_[i + 6] = data_->_des_jacc_cmd[i];
  }

  // Contact Jacobian & Uf & Fr_ieq

  _BuildContactMtxVect();

  _Build_Equality_Constraint();

  _Build_Inequality_Constraint();

  _OptimizationPreparation(Aeq_, beq_, Cieq_, dieq_);

  // Timer timer;
  T f = solve_quadprog(G, g0, CE, ce0, CI, ci0, z);

  // static bool initial(true);
  // if(initial){
  for (size_t i(0); i < dim_opt_; ++i) {
    for (size_t j(0); j < dim_opt_; ++j) {
      g0[i] = G[i][j] * z[j];
    }
  }
  // initial = false;
  //}
  // static int count_ini(0);
  //++count_ini;
  // if(count_ini%100 ==0) std::cout<< "wblc optimization computation: " <<
  // timer.getMs()<<std::endl;

  (void)f;
  _GetSolution(cmd);
  // std::cout << "f: " << f << std::endl;
  // std::cout << "z: " << z << std::endl;
  // char name[1];
  // name[0] = 'x';
  // print_vector(name, z);
  // print_matrix("G", G);
  // name[0] = 'g';
  // print_vector(name, g0);

  // print_matrix("CE", CE);
  // name[0] = 'e';
  // print_vector(name, ce0);

  // print_matrix("CI", CI);
  // name[0] = 'i';
  // print_vector(name, ci0);

  // std::cout << "x: " << z << std::endl;
  // std::cout << "cmd: "<<cmd<<std::endl;

  // if(f > 1.e8){
  // std::cout << "f: " << f << std::endl;
  // std::cout << "x: " << z << std::endl;
  // std::cout << "cmd: "<<cmd<<std::endl;

  // printf("G:\n");
  // std::cout<<G<<std::endl;
  // printf("g0:\n");
  // std::cout<<g0<<std::endl;

  // printf("CE:\n");
  // std::cout<<CE<<std::endl;
  // printf("ce0:\n");
  // std::cout<<ce0<<std::endl;

  // printf("CI:\n");
  // std::cout<<CI<<std::endl;
  // printf("ci0:\n");
  // std::cout<<ci0<<std::endl;
  // exit(0);
  //}
}

template <typename T>
void WBLC<T>::_Build_Inequality_Constraint() {
  size_t row_idx(0);

#if (_WBLC_TORQUE_LIMIT_)
  Cieq_.block(row_idx, WB::num_qdot_, Uf_.rows(), dim_rf_) = Uf_;
  dieq_.head(Uf_.rows()) = Fr_ieq_ - Uf_ * Fr_des_;
  row_idx += Uf_.rows();

  Cieq_.block(row_idx, 0, WB::num_act_joint_, WB::num_qdot_) = WB::Sa_ * WB::A_;
  Cieq_.block(row_idx, WB::num_qdot_, WB::num_act_joint_, dim_rf_) =
      -WB::Sa_ * Jc_.transpose();
  dieq_.segment(row_idx, WB::num_act_joint_) =
      data_->tau_min_ - WB::Sa_ * (WB::cori_ + WB::grav_ + WB::A_ * qddot_);
  row_idx += WB::num_act_joint_;

  Cieq_.block(row_idx, 0, WB::num_act_joint_, WB::num_qdot_) =
      -WB::Sa_ * WB::A_;
  Cieq_.block(row_idx, WB::num_qdot_, WB::num_act_joint_, dim_rf_) =
      WB::Sa_ * Jc_.transpose();
  dieq_.segment(row_idx, WB::num_act_joint_) =
      -data_->tau_max_ + WB::Sa_ * (WB::cori_ + WB::grav_ + WB::A_ * qddot_);
#else
  Cieq_.block(row_idx, 6, Uf_.rows(), dim_rf_) = Uf_;
  dieq_.head(Uf_.rows()) = Fr_ieq_ - Uf_ * Fr_des_;
  row_idx += Uf_.rows();
#endif
  // pretty_print(Cieq_, std::cout, "C ieq");
  // pretty_print(dieq_, std::cout, "d ieq");
}

template <typename T>
void WBLC<T>::_Build_Equality_Constraint() {
#if (_WBLC_TORQUE_LIMIT_)
  // passive joint
  Aeq_.block(0, 0, 6, WB::num_qdot_) = WB::Sv_ * WB::A_;
  Aeq_.block(0, WB::num_qdot_, 6, dim_rf_) = -WB::Sv_ * Jc_.transpose();
  beq_.head(6) = -WB::Sv_ * (WB::A_ * qddot_ + WB::cori_ + WB::grav_ -
                             Jc_.transpose() * Fr_des_);

  // contact xddot
  Aeq_.block(6, 0, dim_rf_, WB::num_qdot_) = Jc_;
  Aeq_.bottomRightCorner(dim_rf_, dim_rf_) =
      -DMat<T>::Identity(dim_rf_, dim_rf_);
  beq_.tail(dim_rf_) = -Jc_ * qddot_ - JcDotQdot_;
#else
  // passive joint
  Aeq_.block(0, 0, 6, 6) = WB::A_.block(0, 0, 6, 6);
  Aeq_.block(0, 6, 6, dim_rf_) = -WB::Sv_ * Jc_.transpose();
  beq_.head(6) = -WB::Sv_ * (WB::A_ * qddot_ + WB::cori_ + WB::grav_ -
                             Jc_.transpose() * Fr_des_);

  // contact xddot
  Aeq_.block(6, 0, dim_rf_, 6) = Jc_.block(0, 0, dim_rf_, 6);
  Aeq_.bottomRightCorner(dim_rf_, dim_rf_) =
      -DMat<T>::Identity(dim_rf_, dim_rf_);
  beq_.tail(dim_rf_) = -Jc_ * qddot_ - JcDotQdot_;
#endif

  // pretty_print(WB::Sv_, std::cout, "Sv");
  // pretty_print(WB::A_, std::cout, "A");
  // pretty_print(WB::cori_, std::cout, "cori");
  // pretty_print(WB::grav_, std::cout, "grav");

  // pretty_print(qddot_, std::cout, "qddot");
  // pretty_print(Jc_, std::cout, "Jc");
  // pretty_print(JcDotQdot_, std::cout, "JcDotQdot");

  // pretty_print(Aeq_, std::cout, "Aeq");
  // pretty_print(beq_, std::cout, "beq");
}

template <typename T>
void WBLC<T>::_BuildContactMtxVect() {
  DMat<T> Jc_i, Uf_i;
  DVec<T> Fr_ieq_i, JcDotQdot_i;

  (_contact_list)[0]->getContactJacobian(Jc_i);
  (_contact_list)[0]->getJcDotQdot(JcDotQdot_i);
  (_contact_list)[0]->getRFConstraintMtx(Uf_i);
  (_contact_list)[0]->getRFConstraintVec(Fr_ieq_i);

  size_t dim_rf = Jc_i.rows();
  size_t num_rows_Uf = Uf_i.rows();

  Jc_.block(0, 0, dim_rf, WB::num_qdot_) = Jc_i;
  JcDotQdot_.head(dim_rf) = JcDotQdot_i;

  Uf_.block(0, 0, num_rows_Uf, dim_rf) = Uf_i;
  Fr_ieq_.head(num_rows_Uf) = Fr_ieq_i;

  Fr_des_.head(dim_rf) = _contact_list[0]->getRFDesired();

  for (size_t i(1); i < (_contact_list).size(); ++i) {
    (_contact_list)[i]->getContactJacobian(Jc_i);
    (_contact_list)[i]->getJcDotQdot(JcDotQdot_i);
    (_contact_list)[i]->getRFConstraintMtx(Uf_i);
    (_contact_list)[i]->getRFConstraintVec(Fr_ieq_i);

    Jc_.block(dim_rf, 0, Jc_i.rows(), WB::num_qdot_) = Jc_i;

    JcDotQdot_.segment(dim_rf, Jc_i.rows()) = JcDotQdot_i;

    Uf_.block(num_rows_Uf, dim_rf, Uf_i.rows(), Uf_i.cols()) = Uf_i;

    Fr_ieq_.segment(num_rows_Uf, Uf_i.rows()) = Fr_ieq_i;

    Fr_des_.segment(dim_rf, Jc_i.rows()) = _contact_list[i]->getRFDesired();

    dim_rf += Jc_i.rows();
    num_rows_Uf += Uf_i.rows();
  }
  // pretty_print(Fr_ieq_, std::cout, "Fr_ieq");
  // pretty_print(Fr_des_, std::cout, "Fr_des");
  // pretty_print(Jc_, std::cout, "Jc");
  // pretty_print(Uf_, std::cout, "Uf");
  // pretty_print(JcDotQdot_, std::cout, "JcDotQdot");
  // printf("size JcDotQdot : %ld, %ld\n", JcDotQdot_.rows(),
  // JcDotQdot_.cols());
}

template <typename T>
void WBLC<T>::_OptimizationPreparation(const DMat<T>& Aeq, const DVec<T>& beq,
                                       const DMat<T>& Cieq,
                                       const DVec<T>& dieq) {
  // Set Cost

#if (_WBLC_TORQUE_LIMIT_)
  for (size_t i(0); i < WB::num_qdot_; ++i) {
    G[i][i] = data_->W_qddot_[i];
  }
  size_t idx_offset = WB::num_qdot_;
#else
  for (size_t i(0); i < 6; ++i) {
    G[i][i] = data_->W_qddot_[i];
  }
  size_t idx_offset = 6;
#endif
  for (size_t i(0); i < dim_rf_; ++i) {
    G[i + idx_offset][i + idx_offset] = data_->W_rf_[i];
  }
  idx_offset += dim_rf_;
  for (size_t i(0); i < dim_rf_; ++i) {
    G[i + idx_offset][i + idx_offset] = data_->W_xddot_[i];
  }

  for (size_t i(0); i < dim_eq_cstr_; ++i) {
    for (size_t j(0); j < dim_opt_; ++j) {
      CE[j][i] = Aeq(i, j);
    }
    ce0[i] = -beq[i];
  }

  for (size_t i(0); i < dim_ieq_cstr_; ++i) {
    for (size_t j(0); j < dim_opt_; ++j) {
      CI[j][i] = Cieq(i, j);
    }
    ci0[i] = -dieq[i];
  }
  // pretty_print(data_->W_rf_, std::cout, "W rf");
  // pretty_print(data_->W_xddot_, std::cout, "W xddot");
  // printf("G:\n");
  // std::cout<<G<<std::endl;
  // printf("g0:\n");
  // std::cout<<g0<<std::endl;

  // printf("CE:\n");
  // std::cout<<CE<<std::endl;
  // printf("ce0:\n");
  // std::cout<<ce0<<std::endl;

  // std::cout<<CI<<std::endl;
  // print_matrix("CI", CI);
  // char name[10];
  // name[0] = 'c';
  // print_vector(name, ci0);
  // std::cout<<ci0<<std::endl;
}

template <typename T>
void WBLC<T>::_GetSolution(DVec<T>& cmd) {
  data_->Fr_ = DVec<T>(dim_rf_);
#if (_WBLC_TORQUE_LIMIT_)
  for (size_t i(0); i < WB::num_qdot_; ++i) qddot_[i] += z[i];
  for (size_t i(0); i < dim_rf_; ++i)
    data_->Fr_[i] = Fr_des_[i] + z[i + WB::num_qdot_];
#else
  for (size_t i(0); i < 6; ++i) qddot_[i] += z[i];
  for (size_t i(0); i < dim_rf_; ++i) data_->Fr_[i] = Fr_des_[i] + z[i + 6];
#endif
  DVec<T> tau =
      WB::A_ * qddot_ + WB::cori_ + WB::grav_ - Jc_.transpose() * data_->Fr_;

  data_->qddot_ = qddot_;
  cmd = tau.tail(WB::num_act_joint_);

  // pretty_print(qddot_, std::cout, "qddot_");
  // pretty_print(delta_qddot, std::cout, "delta_qddot");
  // pretty_print(Fr_des_, std::cout, "Fr des");
  // pretty_print(data_->Fr_, std::cout, "Fr");
  // std::cout<<"z: "<<z<<std::endl;
  // pretty_print(tau, std::cout, "total tau");
  // DVec<T> x_check = Jc_ * (qddot_ + delta_qddot)  + JcDotQdot_;
  // pretty_print(x_check, std::cout, "x check");
}

template class WBLC<double>;
template class WBLC<float>;
