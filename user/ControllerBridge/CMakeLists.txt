include_directories(${PROJECT_SOURCE_DIR})
include_directories("../../robot/include")
include_directories("../../common/include/")
include_directories("../../lcm-types/cpp")

add_executable(ctrl_bridge main.cpp ControllerBridge.cpp)
target_link_libraries(ctrl_bridge robot biomimetics dynacore_param_handler qpOASES)
