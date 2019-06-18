
    if(WIN32)
    set(HPMPC_FLAGS -DOS_WINDOWS)
    elseif(APPLE)
     set(HPMPC_FLAGS -DOS_MAC)
    else()
     set(HPMPC_FLAGS -DOS_LINUX)
    endif()
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../blasfeo/src/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include/foo/bar)
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
add_library(casadi_hpmpc_X64_AVX_GENERIC SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/avx/kernel_dgemm_avx_lib4.c
src/kernel/avx/kernel_dtrmm_avx_lib4.c
src/kernel/avx/kernel_dtrsm_avx_lib4.c
src/kernel/avx/kernel_dsyrk_avx_lib4.c
src/kernel/avx/kernel_dpotrf_avx_lib4.c
src/kernel/avx/kernel_dgemv_avx_lib4.c
src/kernel/avx/kernel_dtrmv_avx_lib4.c
src/kernel/avx/kernel_dtrsv_avx_lib4.c
src/kernel/avx/kernel_dsymv_avx_lib4.c
src/kernel/avx/kernel_dtran_avx_lib4.c
src/kernel/avx/kernel_dttmm_avx_lib4.c
src/kernel/avx/kernel_dtrinv_avx_lib4.c
src/kernel/avx/kernel_dcopy_avx_lib4.c
src/kernel/avx/kernel_dgetrf_avx_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/avx/d_aux_ip_hard_lib4.c
src/mpc_solvers/avx/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX_GENERIC PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_GENERIC -DTARGET_X64_AVX -O2 -fPIC -DBLASFEO -m64 -mavx -DTARGET_X64_AVX ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_X64_AVX_GENERIC
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)
    
target_link_libraries(casadi_hpmpc_X64_AVX_GENERIC casadi_blasfeo_GENERIC m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
add_library(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/avx2/kernel_dgemm_avx2_lib4.c
src/kernel/avx2/kernel_dtrmm_avx2_lib4.c
src/kernel/avx2/kernel_dtrsm_avx2_lib4.c
src/kernel/avx2/kernel_dsyrk_avx2_lib4.c
src/kernel/avx2/kernel_dpotrf_avx2_lib4.c
src/kernel/avx2/kernel_dgemv_avx2_lib4.c
src/kernel/avx2/kernel_dtrmv_avx2_lib4.c
src/kernel/avx2/kernel_dtrsv_avx2_lib4.c
src/kernel/avx2/kernel_dsymv_avx2_lib4.c
src/kernel/avx2/kernel_dtran_avx2_lib4.c
src/kernel/avx2/kernel_dttmm_avx2_lib4.c
src/kernel/avx2/kernel_dtrinv_avx2_lib4.c
src/kernel/avx2/kernel_dgetrf_avx2_lib4.c
src/kernel/avx/kernel_dcopy_avx_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/avx/d_aux_ip_hard_lib4.c
src/mpc_solvers/avx/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_HASWELL -DTARGET_X64_AVX2 -O2 -fPIC -DBLASFEO -m64 -mavx2 -mfma -DTARGET_X64_AVX2 ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)

target_link_libraries(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL casadi_blasfeo_X64_INTEL_HASWELL m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
add_library(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/avx/kernel_dgemm_avx_lib4.c
src/kernel/avx/kernel_dtrmm_avx_lib4.c
src/kernel/avx/kernel_dtrsm_avx_lib4.c
src/kernel/avx/kernel_dsyrk_avx_lib4.c
src/kernel/avx/kernel_dpotrf_avx_lib4.c
src/kernel/avx/kernel_dgemv_avx_lib4.c
src/kernel/avx/kernel_dtrmv_avx_lib4.c
src/kernel/avx/kernel_dtrsv_avx_lib4.c
src/kernel/avx/kernel_dsymv_avx_lib4.c
src/kernel/avx/kernel_dtran_avx_lib4.c
src/kernel/avx/kernel_dttmm_avx_lib4.c
src/kernel/avx/kernel_dtrinv_avx_lib4.c
src/kernel/avx/kernel_dcopy_avx_lib4.c
src/kernel/avx/kernel_dgetrf_avx_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/avx/d_aux_ip_hard_lib4.c
src/mpc_solvers/avx/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_SANDY_BRIDGE -DTARGET_X64_AVX -O2 -fPIC -DBLASFEO -m64 -mavx -DTARGET_X64_AVX ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)

target_link_libraries(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE casadi_blasfeo_X64_INTEL_SANDY_BRIDGE m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
add_library(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/sse3/kernel_dgemm_sse3_lib4.c
src/kernel/sse3/kernel_dtrmm_sse3_lib4.c
src/kernel/sse3/kernel_dtrsm_sse3_lib4.c
src/kernel/sse3/kernel_dsyrk_sse3_lib4.c
src/kernel/sse3/kernel_dpotrf_sse3_lib4.c
src/kernel/c99/kernel_dgemv_c99_lib4.c
src/kernel/c99/kernel_dtrmv_c99_lib4.c
src/kernel/c99/kernel_dtrsv_c99_lib4.c
src/kernel/c99/kernel_dsymv_c99_lib4.c
src/kernel/c99/kernel_dtran_c99_lib4.c
src/kernel/c99/kernel_dttmm_c99_lib4.c
src/kernel/c99/kernel_dtrinv_c99_lib4.c
src/kernel/c99/kernel_dcopy_c99_lib4.c
src/kernel/c99/kernel_dgetrf_c99_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/c99/d_aux_ip_hard_lib4.c
src/mpc_solvers/c99/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_CORE -DTARGET_X64_SSE3 -O2 -fPIC -DBLASFEO -m64 -msse3 -DTARGET_X64_SSE3 ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_X64_SSE3_X64_INTEL_CORE
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)

target_link_libraries(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE casadi_blasfeo_X64_INTEL_CORE m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
add_library(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/sse3/kernel_dgemm_sse3_lib4.c
src/kernel/sse3/kernel_dtrmm_sse3_lib4.c
src/kernel/sse3/kernel_dtrsm_sse3_lib4.c
src/kernel/sse3/kernel_dsyrk_sse3_lib4.c
src/kernel/sse3/kernel_dpotrf_sse3_lib4.c
src/kernel/c99/kernel_dgemv_c99_lib4.c
src/kernel/c99/kernel_dtrmv_c99_lib4.c
src/kernel/c99/kernel_dtrsv_c99_lib4.c
src/kernel/c99/kernel_dsymv_c99_lib4.c
src/kernel/c99/kernel_dtran_c99_lib4.c
src/kernel/c99/kernel_dttmm_c99_lib4.c
src/kernel/c99/kernel_dtrinv_c99_lib4.c
src/kernel/c99/kernel_dcopy_c99_lib4.c
src/kernel/c99/kernel_dgetrf_c99_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/c99/d_aux_ip_hard_lib4.c
src/mpc_solvers/c99/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_AMD_BULLDOZER -DTARGET_X64_SSE3 -O2 -fPIC -DBLASFEO -m64 -msse3 -DTARGET_X64_SSE3 ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)

target_link_libraries(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER casadi_blasfeo_X64_AMD_BULLDOZER m)
endif()
add_library(casadi_hpmpc_C99_4X4_GENERIC SHARED
src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/kernel/c99/kernel_dgemm_c99_lib4.c
src/kernel/c99/kernel_dtrmm_c99_lib4.c
src/kernel/c99/kernel_dtrsm_c99_lib4.c
src/kernel/c99/kernel_dsyrk_c99_lib4.c
src/kernel/c99/kernel_dpotrf_c99_lib4.c
src/kernel/c99/kernel_dgemv_c99_lib4.c
src/kernel/c99/kernel_dtrmv_c99_lib4.c
src/kernel/c99/kernel_dtrsv_c99_lib4.c
src/kernel/c99/kernel_dsymv_c99_lib4.c
src/kernel/c99/kernel_dtran_c99_lib4.c
src/kernel/c99/kernel_dttmm_c99_lib4.c
src/kernel/c99/kernel_dtrinv_c99_lib4.c
src/kernel/c99/kernel_dcopy_c99_lib4.c
src/kernel/c99/kernel_dgetrf_c99_lib4.c
src/blas/blas_d_lib4.c
src/lqcp_solvers/d_back_ric_rec.c
src/lqcp_solvers/d_for_schur_rec.c
src/lqcp_solvers/d_res.c
src/lqcp_solvers/d_part_cond.c
src/mpc_solvers/d_ip2_hard.c
src/mpc_solvers/d_res_ip_hard.c
src/mpc_solvers/d_ip2_res_hard.c
src/mpc_solvers/c99/d_aux_ip_hard_lib4.c
src/mpc_solvers/c99/d_res_ip_res_hard.c
src/interfaces/c/c_interface_work_space.c
src/interfaces/c/c_order_interface.c
src/interfaces/c/fortran_order_interface.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_C99_4X4_GENERIC PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_GENERIC -DTARGET_C99_4X4 -O2 -fPIC -DBLASFEO -march=native -DTARGET_C99_4X4 ${HPMPC_FLAGS}")
install(TARGETS casadi_hpmpc_C99_4X4_GENERIC
  LIBRARY DESTINATION ${LIB_PREFIX}
  ARCHIVE DESTINATION ${LIB_PREFIX}
  RUNTIME DESTINATION ${BIN_PREFIX}
)

target_link_libraries(casadi_hpmpc_C99_4X4_GENERIC casadi_blasfeo_GENERIC m)
