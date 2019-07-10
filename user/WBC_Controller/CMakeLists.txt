include_directories("./")
include_directories("../../robot/include")
include_directories("../../common/include/")
include_directories("../../third-party/")
include_directories("../../third-party/ParamHandler")

include_directories("../../lcm-types/cpp")
include_directories("/usr/local/include/lcm/")   # lcm includes
include_directories(${CMAKE_BINARY_DIR})


add_subdirectory(WBC_States)
add_subdirectory(WBC)

add_executable(wbc_ctrl main.cpp WBC_Controller.cpp)
target_link_libraries(wbc_ctrl robot biomimetics)
target_link_libraries(wbc_ctrl WBC_state)

