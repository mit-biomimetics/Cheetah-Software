cmake_minimum_required(VERSION 2.8.6)

# Linear solver using dense LU factorization in LAPACK
casadi_plugin(Linsol lapacklu
  lapack_lu.hpp
  lapack_lu.cpp
  lapack_lu_meta.cpp)
casadi_plugin_link_libraries(Linsol lapacklu ${LAPACK_LIBRARIES})

# Linear solver using dense QR factorization in LAPACK
casadi_plugin(Linsol lapackqr
  lapack_qr.hpp
  lapack_qr.cpp
  lapack_qr_meta.cpp
  )
casadi_plugin_link_libraries(Linsol lapackqr ${LAPACK_LIBRARIES})
