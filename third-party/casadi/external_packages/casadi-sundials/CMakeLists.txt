# This is a alternative version of the original CMakeLists.txt
# Joel Andersson
# K.U. Leuven 2012, 2016

cmake_minimum_required(VERSION 2.8.6)

project(sundials-2.6.1mod C)
include_directories(.)

set(SUNDIALS_INCLUDES
  nvector/nvector_serial.h

  sundials/sundials_band.h
  sundials/sundials_dense.h
  sundials/sundials_direct.h
  sundials/sundials_iterative.h
#   sundials/sundials_lapack.h
  sundials/sundials_math.h
  sundials/sundials_nvector.h
  sundials/sundials_spbcgs.h
  sundials/sundials_spgmr.h
  sundials/sundials_sptfqmr.h
  sundials/sundials_types.h

  cvodes/cvodes_band.h
  cvodes/cvodes_bandpre.h
  cvodes/cvodes_bbdpre.h
  cvodes/cvodes_dense.h
  cvodes/cvodes_diag.h
  cvodes/cvodes_direct.h
  cvodes/cvodes.h
#   cvodes/cvodes_lapack.h
  cvodes/cvodes_spbcgs.h
  cvodes/cvodes_spgmr.h
  cvodes/cvodes_spils.h
  cvodes/cvodes_sptfqmr.h

  idas/idas_band.h
  idas/idas_bbdpre.h
  idas/idas_dense.h
  idas/idas_direct.h
  idas/idas.h
#   idas/idas_lapack.h
  idas/idas_spbcgs.h
  idas/idas_spgmr.h
  idas/idas_spils.h
  idas/idas_sptfqmr.h

  kinsol/kinsol_band.h
  kinsol/kinsol_bbdpre.h
  kinsol/kinsol_dense.h
  kinsol/kinsol_direct.h
  kinsol/kinsol.h
#   kinsol/kinsol_lapack.h
  kinsol/kinsol_spbcgs.h
  kinsol/kinsol_spgmr.h
  kinsol/kinsol_spils.h
  kinsol/kinsol_sptfqmr.h
)

set(SUNDIALS_SOURCES
  nvector/nvector_serial.c

  sundials/sundials_band.c
  sundials/sundials_dense.c
  sundials/sundials_direct.c
  sundials/sundials_iterative.c
  sundials/sundials_math.c
  sundials/sundials_nvector.c
  sundials/sundials_spbcgs.c
  sundials/sundials_spgmr.c
  sundials/sundials_sptfqmr.c

  cvodes/cvodea.c
  cvodes/cvodea_io.c
  cvodes/cvodes_band.c
  cvodes/cvodes_bandpre.c
  cvodes/cvodes_bandpre_impl.h
  cvodes/cvodes_bbdpre.c
  cvodes/cvodes_bbdpre_impl.h
  cvodes/cvodes.c
  cvodes/cvodes_dense.c
  cvodes/cvodes_diag.c
  cvodes/cvodes_diag_impl.h
  cvodes/cvodes_direct.c
  cvodes/cvodes_direct_impl.h
  cvodes/cvodes_impl.h
  cvodes/cvodes_io.c
#   cvodes/cvodes_lapack.c
  cvodes/cvodes_spbcgs.c
  cvodes/cvodes_spgmr.c
  cvodes/cvodes_spils.c
  cvodes/cvodes_spils_impl.h
  cvodes/cvodes_sptfqmr.c

  idas/idaa.c
  idas/idaa_io.c
  idas/idas_band.c
  idas/idas_bbdpre.c
  idas/idas_bbdpre_impl.h
  idas/idas.c
  idas/idas_dense.c
  idas/idas_direct.c
  idas/idas_direct_impl.h
  idas/idas_ic.c
  idas/idas_impl.h
  idas/idas_io.c
#   idas/idas_lapack.c
  idas/idas_spbcgs.c
  idas/idas_spgmr.c
  idas/idas_spils.c
  idas/idas_spils_impl.h
  idas/idas_sptfqmr.c

  kinsol/kinsol_band.c
  kinsol/kinsol_bbdpre.c
  kinsol/kinsol_bbdpre_impl.h
  kinsol/kinsol.c
  kinsol/kinsol_dense.c
  kinsol/kinsol_direct.c
  kinsol/kinsol_direct_impl.h
  kinsol/kinsol_impl.h
  kinsol/kinsol_io.c
#   kinsol/kinsol_lapack.c
  kinsol/kinsol_spbcgs.c
  kinsol/kinsol_spgmr.c
  kinsol/kinsol_spils.c
  kinsol/kinsol_spils_impl.h
  kinsol/kinsol_sptfqmr.c
)

add_library(casadi_sundials STATIC ${SUNDIALS_INCLUDES} ${SUNDIALS_SOURCES})
