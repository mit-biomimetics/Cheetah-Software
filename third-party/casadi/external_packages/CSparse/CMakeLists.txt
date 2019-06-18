cmake_minimum_required(VERSION 2.8.6)

include_directories(Include)

set(CSPARSE_SRC
  Include/cs.h
  Source/cs_add.c Source/cs_amd.c Source/cs_chol.c Source/cs_cholsol.c Source/cs_counts.c Source/cs_cumsum.c
  Source/cs_droptol.c Source/cs_dropzeros.c Source/cs_dupl.c Source/cs_entry.c
  Source/cs_etree.c Source/cs_fkeep.c Source/cs_gaxpy.c Source/cs_happly.c Source/cs_house.c Source/cs_ipvec.c
  Source/cs_lsolve.c Source/cs_ltsolve.c Source/cs_lu.c Source/cs_lusol.c Source/cs_util.c Source/cs_multiply.c
  Source/cs_permute.c Source/cs_pinv.c Source/cs_post.c Source/cs_pvec.c Source/cs_qr.c Source/cs_qrsol.c
  Source/cs_scatter.c Source/cs_schol.c Source/cs_sqr.c Source/cs_symperm.c Source/cs_tdfs.c Source/cs_malloc.c
  Source/cs_transpose.c Source/cs_compress.c Source/cs_usolve.c Source/cs_utsolve.c Source/cs_scc.c
  Source/cs_maxtrans.c Source/cs_dmperm.c Source/cs_updown.c Source/cs_print.c Source/cs_norm.c Source/cs_load.c
  Source/cs_dfs.c Source/cs_reach.c Source/cs_spsolve.c Source/cs_ereach.c Source/cs_leaf.c Source/cs_randperm.c
)

add_library(casadi_csparse STATIC ${CSPARSE_SRC})
