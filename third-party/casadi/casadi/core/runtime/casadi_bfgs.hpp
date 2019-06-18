// NOLINT(legal/copyright)
// SYMBOL "bfgs"
// BFGS update
template<typename T1>
void casadi_bfgs(const casadi_int* sp_h, T1* h, const T1* dx,
                 const T1* glag, const T1* glag_old, T1* w) {
  // Local variables
  casadi_int nx;
  T1 *yk, *qk, dxBkdx, omega, theta, phi;
  // Dimension
  nx = sp_h[0];
  // Work vectors
  yk = w; w += nx;
  qk = w; w += nx;
  // yk = glag - glag_old
  casadi_copy(glag, nx, yk);
  casadi_axpy(nx, -1., glag_old, yk);
  // qk = H*dx
  casadi_fill(qk, nx, 0.);
  casadi_mv(h, sp_h, dx, qk, 0);
  // Calculating theta
  dxBkdx = casadi_dot(nx, dx, qk);
  // C-REPLACE "if_else" "casadi_if_else"
  omega = if_else(casadi_dot(nx, yk, dx) < 0.2 * casadi_dot(nx, dx, qk),
                  0.8 * dxBkdx / (dxBkdx - casadi_dot(nx, dx, yk)), 1);
  // yk = omega * yk + (1 - omega) * qk;
  casadi_scal(nx, omega, yk);
  casadi_axpy(nx, 1 - omega, qk, yk);
  theta = 1. / casadi_dot(nx, dx, yk);
  phi = 1. / casadi_dot(nx, qk, dx);
  // Update H
  casadi_rank1(h, sp_h, theta, yk, yk);
  casadi_rank1(h, sp_h, -phi, qk, qk);
}

// SYMBOL "bfgs_reset"
// Removes off-diagonal entries
template<typename T1>
void casadi_bfgs_reset(const casadi_int* sp_h, T1* h) {
  casadi_int ncol, c, k;
  const casadi_int *colind, *row;
  ncol = sp_h[1];
  colind = sp_h+2; row = sp_h+ncol+3;
  for (c=0; c<ncol; ++c) {
    for (k=colind[c]; k<colind[c+1]; ++k) {
      if (c!=row[k]) h[k] = 0;
    }
  }
}
