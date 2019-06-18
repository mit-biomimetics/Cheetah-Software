# Add IPOPT link directories
if(WITH_QPOASES AND WITH_IPOPT)
  add_executable(shallow_water shallow_water.cpp)
  target_link_libraries(shallow_water casadi)
endif()

# Writing a multiple shooting code from scratch
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(blocksqp_multiple_shooting blocksqp_multiple_shooting.cpp)
  target_link_libraries(blocksqp_multiple_shooting casadi)

  add_executable(blocksqp_test blocksqp_test.cpp)
  target_link_libraries(blocksqp_test casadi)
endif()
