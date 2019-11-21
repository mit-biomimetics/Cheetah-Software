#include "proj.h"


void project(OSQPWorkspace *work, c_float *z) {
  c_int i, m;

  m = work->data->m;

  for (i = 0; i < m; i++) {
    z[i] = c_min(c_max(z[i],
                       work->data->l[i]), // Between lower
                 work->data->u[i]);       // and upper bounds
  }
}

void project_normalcone(OSQPWorkspace *work, c_float *z, c_float *y) {
  c_int i, m;

  // NB: Use z_prev as temporary vector

  m = work->data->m;

  for (i = 0; i < m; i++) {
    work->z_prev[i] = z[i] + y[i];
    z[i]            = c_min(c_max(work->z_prev[i], work->data->l[i]),
                            work->data->u[i]);
    y[i] = work->z_prev[i] - z[i];
  }
}
