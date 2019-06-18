// NOLINT(legal/copyright)
// SYMBOL "nd_boor_eval"
template<typename T1>
void casadi_nd_boor_eval(T1* ret, casadi_int n_dims, const T1* all_knots, const casadi_int* offset, const casadi_int* all_degree, const casadi_int* strides, const T1* c, casadi_int m, const T1* all_x, const casadi_int* lookup_mode, casadi_int reverse, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int n_iter, k, i, pivot;
  casadi_int *boor_offset, *starts, *index, *coeff_offset;
  T1 *cumprod, *all_boor;

  boor_offset = iw; iw+=n_dims+1;
  starts = iw; iw+=n_dims;
  index = iw; iw+=n_dims;
  coeff_offset = iw;

  cumprod = w; w+= n_dims+1;
  all_boor = w;

  boor_offset[0] = 0;
  cumprod[n_dims] = 1;
  coeff_offset[n_dims] = 0;

  n_iter = 1;
  for (k=0;k<n_dims;++k) {
    T1 *boor;
    const T1* knots;
    T1 x;
    casadi_int degree, n_knots, n_b, L, start;
    boor = all_boor+boor_offset[k];

    degree = all_degree[k];
    knots = all_knots + offset[k];
    n_knots = offset[k+1]-offset[k];
    n_b = n_knots-degree-1;

    x = all_x[k];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[k]);

    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;

    starts[k] = start;

    casadi_fill(boor, 2*degree+1, 0.0);
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        casadi_fill(boor, degree+1, 1.0);
      } else if (x==knots[n_knots-1]) {
        boor[degree] = 1;
      } else if (knots[L+degree]==x) {
        boor[degree-1] = 1;
      } else {
        boor[degree] = 1;
      }
    }
    casadi_de_boor(x, knots+start, 2*degree+2, degree, boor);
    boor+= degree+1;
    n_iter*= degree+1;
    boor_offset[k+1] = boor_offset[k] + degree+1;
  }

  casadi_fill_casadi_int(index, n_dims, 0);

  // Prepare cumulative product
  for (pivot=n_dims-1;pivot>=0;--pivot) {
    cumprod[pivot] = (*(all_boor+boor_offset[pivot]))*cumprod[pivot+1];
    coeff_offset[pivot] = starts[pivot]*strides[pivot]+coeff_offset[pivot+1];
  }

  for (k=0;k<n_iter;++k) {
    casadi_int pivot = 0;
    // accumulate result
    for (i=0;i<m;++i) {
      if (reverse) {
        ret[coeff_offset[0]+i] += c[i]*cumprod[0];
      } else {
        ret[i] += c[coeff_offset[0]+i]*cumprod[0];
      }
    }

    // Increment index
    index[0]++;

    // Handle index overflow
    {
      // increment next index (forward)
      while (index[pivot]==boor_offset[pivot+1]-boor_offset[pivot]) {
        index[pivot] = 0;
        if (pivot==n_dims-1) break;
        index[++pivot]++;
      }

      // update cumulative structures (reverse)
      while (pivot>0) {
        // Compute product
        cumprod[pivot] = (*(all_boor+boor_offset[pivot]+index[pivot]))*cumprod[pivot+1];
        // Compute offset
        coeff_offset[pivot] = (starts[pivot]+index[pivot])*strides[pivot]+coeff_offset[pivot+1];
        pivot--;
      }
    }

    // Compute product
    cumprod[0] = (*(all_boor+index[0]))*cumprod[1];

    // Compute offset
    coeff_offset[0] = (starts[0]+index[0])*m+coeff_offset[1];

  }
}
