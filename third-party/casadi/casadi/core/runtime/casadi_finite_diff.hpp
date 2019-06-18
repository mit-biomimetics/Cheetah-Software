// NOLINT(legal/copyright)
// SYMBOL "finite_diff_mem"
template<typename T1>
struct casadi_finite_diff_mem {
  // Input precision
  T1 reltol;
  // Output precision
  T1 abstol;
  // Smoothness parameter
  // Smaller epsilon: More discontinuity rejecting
  // Larger epsilon: More accurate (higher order) if smooth
  T1 smoothing;
};

// C-REPLACE "casadi_finite_diff_mem<T1>" "struct casadi_finite_diff_mem"
// C-REPLACE "std::numeric_limits<T1>::quiet_NaN()" "NAN"
// C-REPLACE "fmax" "casadi_fmax"

// SYMBOL "forward_diff"
template<typename T1>
T1 casadi_forward_diff(T1** yk, T1* y0, T1* J,
                       T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  casadi_int i;
  for (i=0; i<n_y; ++i) {
    J[i] = (yk[0][i]-y0[i])/h;
  }
  return -1;
}

// SYMBOL "central_diff"
template<typename T1>
T1 casadi_central_diff(T1** yk, T1* y0, T1* J,
                       T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  // Return value
  T1 u;
  // Stencil
  T1 yf, yc, yb;
  // Local variables
  T1 err_trunc, err_round;
  casadi_int i;
  // Set u and stencils to zero (also supresses warnings)
  yf = yc = yb = u = 0;
  for (i=0; i<n_y; ++i) {
    // Copy to local variables, return -1 if invalid entry
    if (!isfinite((yf=yk[1][i])) || !isfinite((yc=y0[i])) || !isfinite((yb=yk[0][i]))) {
      J[i] = std::numeric_limits<T1>::quiet_NaN();
      u = -1;
      continue;
    }
    // Central difference approximation
    J[i] = (yf - yb)/(2*h);
    // Truncation error
    err_trunc = yf - 2*yc + yb;
    // Roundoff error
    err_round = m->reltol/h*fmax(fabs(yf-yc), fabs(yc-yb)) + m->abstol;
    // Update error estimate
    if (u>=0) u = fmax(u, fabs(err_trunc/err_round));
  }
  return u;
}

// SYMBOL "smoothing_diff"
template<typename T1>
T1 casadi_smoothing_diff(T1** yk, T1* y0, T1* J,
                         T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  // Return value
  T1 u;
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 Jk, wk, sw, ui, err_trunc, err_round, sm;
  casadi_int i, k;
  // Set u and stencils to zero (also supresses warnings)
  yf = yc = yb = u = 0;
  for (i=0; i<n_y; ++i) {
    // Reset derivative estimate, sum of weights, error estimate
    J[i] = sw = ui = 0;
    // For backward shifted, central and forward shifted
    for (k=0; k<3; ++k) {
      // Calculate shifted finite difference approximation
      if (k==0) {
        // Backward shifted
        // 7.10 in Conte & Carl de Boor: Elementary Numerical Analysis (1972)
        // and 25.3.4 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
        if (!isfinite((yc=yk[0][i]))) continue;
        if (!isfinite((yb=yk[1][i]))) continue;
        yf = y0[i];
        Jk = 3*yf - 4*yc + yb;
        wk = 1;
      } else if (k==1) {
        // Central
        // We give this the "nominal weight" 4 since if all weights are equal,
        // this would amount to a five point formula for the derivative
        // (yb2 - 8*yb + 8*yf - y_f2)/(12*h)
        // cf. 25.3.6 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
        if (!isfinite((yf=yk[2][i]))) continue;
        if (!isfinite((yb=yc))) continue;
        yc = y0[i];
        Jk = yf-yb;
        wk = 4;
      } else {
        // Forward shifted
        if (!isfinite((yc=yf))) continue;
        if (!isfinite((yf=yk[3][i]))) continue;
        yb = y0[i];
        Jk = -3*yb + 4*yc - yf;
        wk = 1;
      }
      // Truncation error
      err_trunc = yf - 2*yc + yb;
      // Roundoff error
      err_round = m->reltol/h*fmax(fabs(yf-yc), fabs(yc-yb)) + m->abstol;
      // We use the second order derivative as a smoothness measure
      sm = err_trunc/(h*h);
      // Modify the weight according to smoothness
      wk /= sm*sm + m->smoothing;
      sw += wk;
      // Added weighted contribution to weight and error
      J[i] += wk * Jk;
      ui += wk * fabs(err_trunc/err_round);
    }
    // If sw is 0, no stencil worked
    if (sw==0) {
      // Set component to 0, return -1
      J[i] = std::numeric_limits<T1>::quiet_NaN();
      u = -1;
    } else {
      // Finalize estimate using the sum of weights and the step length
      J[i] /= 2*h*sw;
      if (u>=0) u = fmax(u, ui/sw);
    }
  }
  return u;
}
