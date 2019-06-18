FILE(GLOB_RECURSE headers *.hpp)
FILE(GLOB_RECURSE sources *.cpp)

set(STATE_NAME WBC)

include_directories("./../")
include_directories("../../lcm-types/cpp")
include_directories("/usr/local/include/lcm/")   # lcm includes


add_library (${STATE_NAME}_state SHARED ${headers} ${sources} )
target_link_libraries (${STATE_NAME}_state biomimetics lcm)
target_link_libraries (${STATE_NAME}_state WBIC WBLC dynacore_param_handler)
