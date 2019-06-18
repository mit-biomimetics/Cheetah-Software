cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Nlpsol ampl
  ampl_interface.hpp
  ampl_interface.cpp
  ampl_interface_meta.cpp)
