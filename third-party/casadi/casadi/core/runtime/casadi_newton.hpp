// NOLINT(legal/copyright)
// SYMBOL "newton_mem"

template<typename T1>
struct casadi_newton_mem {
  casadi_int n;
  T1 abstol;
  T1 abstol_step;
  T1* x;
  T1* g;
  T1* jac_g_x;

  const casadi_int* sp_a;
  const casadi_int* sp_v;
  const casadi_int* sp_r;
  const casadi_int* prinv;
  const casadi_int* pc;

  T1* lin_w;
  T1* lin_v;
  T1* lin_r;
  T1* lin_beta;
};

// C-REPLACE "casadi_newton_mem<T1>" "struct casadi_newton_mem"
// SYMBOL "newton"
template<typename T1>
int casadi_newton(const casadi_newton_mem<T1>* m) {
    // Check tolerance on residual
    if (m->abstol>0 && casadi_norm_inf(m->n, m->g) <= m->abstol) return 1;

    // Factorize J
    casadi_qr(m->sp_a, m->jac_g_x, m->lin_w,
              m->sp_v,  m->lin_v, m->sp_r, m->lin_r, m->lin_beta,
              m->prinv, m->pc);
    // Solve J^(-1) g
    casadi_qr_solve(m->g, 1, 0, m->sp_v, m->lin_v, m->sp_r, m->lin_r, m->lin_beta,
                    m->prinv, m->pc, m->lin_w);

    // Update Xk+1 = Xk - J^(-1) g
    casadi_axpy(m->n, -1., m->g, m->x);

    // Check tolerance on step
    if (m->abstol_step>0 && casadi_norm_inf(m->n, m->g) <= m->abstol_step) return 2;

    // We will need another newton step
    return 0;
}
