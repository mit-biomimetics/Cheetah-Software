cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Dple slicot
  slicot_dple.hpp
  slicot_dple.cpp
  slicot_layer.cpp
  slicot_layer.hpp
  slicot_la.hpp
  slicot_dple_meta.cpp)

casadi_plugin(Expm slicot
  slicot_expm.hpp
  slicot_expm.cpp
  slicot_layer.cpp
  slicot_layer.hpp
  slicot_la.hpp
  slicot_expm_meta.cpp)

casadi_plugin_link_libraries(Dple slicot ${SLICOT_LIBRARIES} ${LAPACK_LIBRARIES})
casadi_plugin_link_libraries(Expm slicot ${SLICOT_LIBRARIES} ${LAPACK_LIBRARIES})
