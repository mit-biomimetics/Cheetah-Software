cmake_minimum_required(VERSION 2.8.6)
include_directories(${KNITRO_INCLUDE_DIR})

casadi_plugin(Nlpsol knitro
  knitro_interface.hpp
  knitro_interface.cpp
  knitro_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol knitro ${KNITRO_LIBRARIES})
