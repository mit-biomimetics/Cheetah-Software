// NOLINT(legal/copyright)
// SYMBOL "project"
template<typename T1>
void casadi_project(const T1* x, const casadi_int* sp_x, T1* y, const casadi_int* sp_y, T1* w) {
  casadi_int ncol_x, ncol_y, i, el;
  const casadi_int *colind_x, *row_x, *colind_y, *row_y;
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;
  ncol_y = sp_y[1];
  colind_y = sp_y+2; row_y = sp_y + 2 + ncol_y+1;
  // Loop over columns of x and y
  for (i=0; i<ncol_x; ++i) {
    // Zero out requested entries in y
    for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;
    // Set x entries
    for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];
    // Retrieve requested entries in y
    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] = w[row_y[el]];
  }
}
