cmake_minimum_required(VERSION 2.8.6)
include_directories(${CPLEX_INCLUDE_DIR})

casadi_plugin(Conic cplex
  cplex_interface.hpp
  cplex_interface.cpp
  cplex_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cplex ${CPLEX_LIBRARIES})

if(CPLEX_DEFINITIONS)
  add_definitions(${CPLEX_DEFINITIONS})
endif()
