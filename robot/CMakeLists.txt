cmake_minimum_required(VERSION 3.5)
project(robot)


include_directories("include/")
include_directories("./")
include_directories("../common/include/")
include_directories("../third-party")
include_directories("../third-party/ParamHandler")
include_directories("../third-party/vectornav/include")
include_directories("../third-party/vectornav/include/vn")
include_directories(SYSTEM "../third-party/qpOASES/include/")

include_directories("../lcm-types/cpp")
include_directories("/usr/local/include/lcm/")   # lcm includes

include_directories(${CMAKE_BINARY_DIR})
file(GLOB sources "src/*.cpp" "src/rt/*.c"  "src/rt/*.cpp" "FSM_States/*.h" "FSM_States/*.cpp" "Controllers/BalanceController/*.cpp" "Controllers/BalanceController/*.h" "Controllers/BalanceController/*.hpp")

add_executable(robot ${sources})
target_link_libraries(robot biomimetics pthread lcm rt inih libvnc)

add_subdirectory(WBC_States)
add_subdirectory(WBC)

target_link_libraries(robot dynacore_param_handler WBC_state qpOASES)

