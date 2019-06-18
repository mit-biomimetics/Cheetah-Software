cmake_minimum_required(VERSION 2.8.6)
include_directories(${CSPARSE_INCLUDE_DIR})

casadi_plugin(Linsol csparse
  csparse_interface.hpp
  csparse_interface.cpp
  csparse_interface_meta.cpp)
casadi_plugin_link_libraries(Linsol csparse ${CSPARSE_LIBRARIES})

casadi_plugin(Linsol csparsecholesky
  csparse_cholesky_interface.hpp
  csparse_cholesky_interface.cpp
  csparse_cholesky_interface_meta.cpp
  )
casadi_plugin_link_libraries(Linsol csparsecholesky ${CSPARSE_LIBRARIES})
