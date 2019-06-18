// NOLINT(legal/copyright)
// SYMBOL "flip"
inline
casadi_int casadi_flip(casadi_int* corner, casadi_int ndim) {
  casadi_int i;
  for (i=0; i<ndim; ++i) {
    if (corner[i]) {
      corner[i]=0;
    } else {
      corner[i]=1;
      return 1;
    }
  }
  return 0;
}
