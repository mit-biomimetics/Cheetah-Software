cmake_minimum_required(VERSION 2.8.6)
include_directories(${SNOPT_INTERFACE_INCLUDE_DIR})

casadi_plugin(Nlpsol snopt
  snopt_interface.hpp
  snopt_interface.cpp
  snopt_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol snopt ${SNOPT_INTERFACE_LIBRARIES} ${SNOPT_LIBRARIES})
