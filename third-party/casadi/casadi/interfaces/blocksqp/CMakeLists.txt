cmake_minimum_required(VERSION 2.8.6)
include_directories(${QPOASES_INCLUDE_DIR})

casadi_plugin(Nlpsol blocksqp
  blocksqp.hpp
  blocksqp.cpp
  blocksqp_meta.cpp)

casadi_plugin_link_libraries(Nlpsol blocksqp ${QPOASES_LIBRARIES})
casadi_plugin_link_libraries(Nlpsol blocksqp casadi_conic_qpoases)
