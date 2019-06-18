// NOLINT(legal/copyright)
// SYMBOL "interpn_interpolate"
template<typename T1>
void casadi_interpn_interpolate(T1* res, casadi_int ndim, const casadi_int* offset, const T1* values, const T1* alpha, const casadi_int* index, const casadi_int* corner, T1* coeff, casadi_int m) { // NOLINT(whitespace/line_length)
  T1 c;
  casadi_int ld, i;
  // Get weight and value for corner
  c=1;
  ld=1; // leading dimension
  for (i=0; i<ndim; ++i) {
    if (coeff) *coeff++ = c;
    if (corner[i]) {
      c *= alpha[i];
    } else {
      c *= 1-alpha[i];
    }
    values += (index[i]+corner[i])*ld*m;
    ld *= offset[i+1]-offset[i];
  }
  if (coeff) {
    for (i=0;i<m;++i) res[i] += values[i];
  } else {
    for (i=0;i<m;++i) res[i] += c*values[i];
  }
}
