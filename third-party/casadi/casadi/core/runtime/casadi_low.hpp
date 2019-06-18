// NOLINT(legal/copyright)
// SYMBOL "low"
template<typename T1>
casadi_int casadi_low(T1 x, const double* grid, casadi_int ng, casadi_int lookup_mode) {
  switch (lookup_mode) {
    case 1: // exact
      {
        double g0, dg;
        casadi_int ret;
        g0 = grid[0];
        dg = grid[ng-1]-g0;
        ret = (casadi_int) ((x-g0)*(ng-1)/dg); // NOLINT(readability/casting)
        if (ret<0) ret=0;
        if (ret>ng-2) ret=ng-2;
        return ret;
      }
    case 2: // binary
      {
        casadi_int start, stop, pivot;
        // Quick return
        if (ng<2 || x<grid[1]) return 0;
        if (x>grid[ng-1]) return ng-2;

        start = 0;
        stop  = ng-1;
        while (1) {
          pivot = (stop+start)/2;
          if (x < grid[pivot]) {
            if (pivot==stop) return pivot;
            stop = pivot;
          } else {
            if (pivot==start) return pivot;
            start = pivot;
          }
        }
      }
    default: // linear
      {
        casadi_int i;
        for (i=0; i<ng-2; ++i) {
          if (x < grid[i+1]) break;
        }
        return i;
      }
  }
}
