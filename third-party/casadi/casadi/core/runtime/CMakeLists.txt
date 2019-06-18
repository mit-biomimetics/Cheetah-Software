cmake_minimum_required(VERSION 2.8.6)

# C runtime sources
set(RUNTIME_SRC
  casadi_axpy.hpp
  casadi_bilin.hpp
  casadi_copy.hpp
  casadi_de_boor.hpp
  casadi_densify.hpp
  casadi_dot.hpp
  casadi_fill.hpp
  casadi_flip.hpp
  casadi_getu.hpp
  casadi_iamax.hpp
  casadi_interpn.hpp
  casadi_interpn_grad.hpp
  casadi_interpn_interpolate.hpp
  casadi_interpn_weights.hpp
  casadi_low.hpp
  casadi_max_viol.hpp
  casadi_minmax.hpp
  casadi_mtimes.hpp
  casadi_mv.hpp
  casadi_mv_dense.hpp
  casadi_nd_boor_eval.hpp
  casadi_norm_1.hpp
  casadi_norm_2.hpp
  casadi_norm_inf.hpp
  casadi_norm_inf_mul.hpp
  casadi_polyval.hpp
  casadi_project.hpp
  casadi_rank1.hpp
  casadi_scal.hpp
  casadi_sparsify.hpp
  casadi_sum_viol.hpp
  casadi_swap.hpp
  casadi_trans.hpp
  casadi_finite_diff.hpp
  casadi_ldl.hpp
  casadi_qr.hpp
  casadi_qp.hpp
  casadi_bfgs.hpp
  casadi_regularize.hpp
  casadi_newton.hpp
)
set(CASADI_RUNTIME_SRC "${RUNTIME_SRC}" PARENT_SCOPE)

# Stringify C runtime
CASADI_STRINGIFY("${CMAKE_CURRENT_BINARY_DIR}/casadi_runtime_str.h"
  ${RUNTIME_SRC}
  casadi_to_mex.hpp
  casadi_from_mex.hpp
)

set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Wno-conversion
-Wno-error")

install(FILES casadi_runtime.hpp shared.hpp ${RUNTIME_SRC}
  DESTINATION ${INCLUDE_PREFIX}/casadi/core/runtime
)
