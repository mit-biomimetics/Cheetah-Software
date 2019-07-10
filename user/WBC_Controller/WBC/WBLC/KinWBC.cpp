#include "KinWBC.hpp"
#include <Utilities/Utilities_print.h>
#include <Utilities/pseudoInverse.h>

template <typename T>
KinWBC<T>::KinWBC(size_t num_qdot)
    : threshold_(0.001), num_qdot_(num_qdot), num_act_joint_(num_qdot - 6) {
  I_mtx = DMat<T>::Identity(num_qdot_, num_qdot_);
}

template <typename T>
bool KinWBC<T>::FindConfiguration(
    const DVec<T>& curr_config, const std::vector<Task<T>*>& task_list,
    const std::vector<ContactSpec<T>*>& contact_list, DVec<T>& jpos_cmd,
    DVec<T>& jvel_cmd, DVec<T>& jacc_cmd) {
  // Contact Jacobian Setup
  DMat<T> Jc, Jc_i;
  contact_list[0]->getContactJacobian(Jc);
  size_t num_rows = Jc.rows();

  for (size_t i(1); i < contact_list.size(); ++i) {
    contact_list[i]->getContactJacobian(Jc_i);
    size_t num_new_rows = Jc_i.rows();
    Jc.conservativeResize(num_rows + num_new_rows, num_qdot_);
    Jc.block(num_rows, 0, num_new_rows, num_qdot_) = Jc_i;
    num_rows += num_new_rows;
  }

  // Projection Matrix
  DMat<T> Nc;
  _BuildProjectionMatrix(Jc, Nc);

  DVec<T> delta_q, qdot, qddot, JtDotQdot;
  DMat<T> Jt, JtPre, JtPre_pinv, N_nx, N_pre;

  // First Task
  Task<T>* task = task_list[0];
  task->getTaskJacobian(Jt);
  task->getTaskJacobianDotQdot(JtDotQdot);
  JtPre = Jt * Nc;
  _PseudoInverse(JtPre, JtPre_pinv);

  delta_q = JtPre_pinv * (task->getPosError());
  qdot = JtPre_pinv * (task->getDesVel());
  qddot = JtPre_pinv * (task->getDesAcc() - JtDotQdot);

  DVec<T> prev_delta_q = delta_q;
  DVec<T> prev_qdot = qdot;
  DVec<T> prev_qddot = qddot;

  _BuildProjectionMatrix(JtPre, N_nx);
  N_pre = Nc * N_nx;

  // pretty_print(qddot, std::cout, "qddot");
  // pretty_print(JtDotQdot, std::cout, "JtDotQdot");

  // DVec<T> xdot_c1 = Jc * delta_q;
  // pretty_print(xdot_c1, std::cout, "1st contact vel");
  // pretty_print(Jt, std::cout, "1st task Jt");
  // pretty_print(Jc, std::cout, "Jc");
  // pretty_print(Nc, std::cout, "Nc");
  // pretty_print(JtPre, std::cout, "JtNc");
  // pretty_print(JtPre_pinv, std::cout, "JtNc_inv");
  // pretty_print(delta_q, std::cout, "delta q");
  // DMat<T> test = Jt * N_pre;
  // pretty_print(test, std::cout, "Jt1N1");

  for (size_t i(1); i < task_list.size(); ++i) {
    task = task_list[i];

    task->getTaskJacobian(Jt);
    task->getTaskJacobianDotQdot(JtDotQdot);
    JtPre = Jt * N_pre;

    _PseudoInverse(JtPre, JtPre_pinv);
    delta_q =
        prev_delta_q + JtPre_pinv * (task->getPosError() - Jt * prev_delta_q);
    qdot = prev_qdot + JtPre_pinv * (task->getDesVel() - Jt * prev_qdot);
    qddot = prev_qddot +
            JtPre_pinv * (task->getDesAcc() - JtDotQdot - Jt * prev_qddot);

    // pretty_print(Jt, std::cout, "2nd Jt");
    // pretty_print(N_pre, std::cout, "N_pre");
    // pretty_print(JtPre, std::cout, "JtPre");
    // pretty_print(JtPre_pinv, std::cout, "JtPre_inv");
    // pretty_print(delta_q, std::cout, "delta q");

    // For the next task
    _BuildProjectionMatrix(JtPre, N_nx);
    N_pre *= N_nx;
    prev_delta_q = delta_q;
    prev_qdot = qdot;
    prev_qddot = qddot;

    // printf("%lu th\n", i);
    // pretty_print(qddot, std::cout, "qddot");
    // pretty_print(JtDotQdot, std::cout, "JtDotQdot");
  }
  // DVec<T> xdot_c = Jc * delta_q;
  // pretty_print(xdot_c, std::cout, "contact vel");
  for (size_t i(0); i < num_act_joint_; ++i) {
    jpos_cmd[i] = curr_config[i + 6] + delta_q[i + 6];
    jvel_cmd[i] = qdot[i + 6];
    jacc_cmd[i] = qddot[i + 6];
  }
  return true;
}

template <typename T>
void KinWBC<T>::_BuildProjectionMatrix(const DMat<T>& J, DMat<T>& N) {
  DMat<T> J_pinv;
  _PseudoInverse(J, J_pinv);
  N = I_mtx - J_pinv * J;
}

template <typename T>
void KinWBC<T>::_PseudoInverse(const DMat<T> J, DMat<T>& Jinv) {
  pseudoInverse(J, threshold_, Jinv);

  // DMat<T> Lambda_inv = J * Ainv_ * J.transpose();
  // DMat<T> Lambda;
  // pseudoInverse(Lambda_inv, threshold_, Lambda);
  // Jinv = Ainv_ * J.transpose() * Lambda;
}

template class KinWBC<float>;
template class KinWBC<double>;
