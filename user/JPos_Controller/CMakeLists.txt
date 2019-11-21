include_directories(${PROJECT_SOURCE_DIR})
include_directories("../../robot/include")
include_directories("../../common/include/")
include_directories("../../lcm-types/cpp")

add_executable(jpos_ctrl main.cpp JPos_Controller.cpp)
target_link_libraries(jpos_ctrl robot biomimetics dynacore_param_handler qpOASES)
