cmake_minimum_required(VERSION 2.8.6)
include_directories(${HPMPC_INCLUDE_DIR})

casadi_plugin(Conic hpmpc
  hpmpc_interface.hpp
  hpmpc_interface.cpp
  hpmpc_interface_meta.cpp)

casadi_plugin_link_libraries(Conic hpmpc ${HPMPC_LIBRARIES} ${BLASFEO_LIBRARIES})

if(HPMPC_DLOPEN)
  add_definitions(-DHPMPC_DLOPEN)
endif()

if(BLASFEO_DLOPEN)
  add_definitions(-DBLASFEO_DLOPEN)
endif()

casadi_plugin_link_libraries(Conic hpmpc ${HPMPC_LIBRARIES} ${BLASFEO_LIBRARIES})
