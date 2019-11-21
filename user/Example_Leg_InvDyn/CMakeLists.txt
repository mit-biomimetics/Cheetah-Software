include_directories("../../robot/include")
include_directories("../../common/include/")
include_directories("../../lcm-types/cpp")

add_executable(leg_invdyn_ctrl main.cpp Leg_InvDyn_Controller.cpp)
target_link_libraries(leg_invdyn_ctrl robot)
