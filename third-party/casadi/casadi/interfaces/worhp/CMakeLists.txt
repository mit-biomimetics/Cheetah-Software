cmake_minimum_required(VERSION 2.8.6)

include_directories(../)
include_directories(${WORHP_INCLUDE_DIR})

casadi_plugin(Nlpsol worhp
  worhp_interface.hpp
  worhp_interface.cpp
  worhp_interface_meta.cpp)
casadi_plugin_link_libraries(Nlpsol worhp ${WORHP_LIBRARIES})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_nlpsol_worhp PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
endif()
