// NOLINT(legal/copyright)
// SYMBOL "copy"
template<typename T1>
void casadi_copy(const T1* x, casadi_int n, T1* y) {
  casadi_int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}
