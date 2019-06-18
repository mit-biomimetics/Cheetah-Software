cmake_minimum_required(VERSION 2.8.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

# SUNDIALS commons
casadi_library(casadi_sundials_common sundials_interface.hpp sundials_interface.cpp)
target_link_libraries(casadi_sundials_common ${SUNDIALS_LIBRARIES})

# CVODES
casadi_plugin(Integrator cvodes
  cvodes_interface.hpp    cvodes_interface.cpp
  cvodes_interface_meta.cpp
  )
casadi_plugin_link_libraries(Integrator cvodes casadi_sundials_common)

# IDAS
casadi_plugin(Integrator idas
  idas_interface.hpp      idas_interface.cpp
  idas_interface_meta.cpp
  )
casadi_plugin_link_libraries(Integrator idas casadi_sundials_common)

# KINSOL
casadi_plugin(Rootfinder kinsol
   kinsol_interface.hpp    kinsol_interface.cpp
   kinsol_interface_meta.cpp
   )
casadi_plugin_link_libraries(Rootfinder kinsol  ${SUNDIALS_LIBRARIES})
