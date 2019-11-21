FILE(GLOB_RECURSE headers *.hpp)
FILE(GLOB_RECURSE sources *.cpp)

add_library (VisionMPC SHARED ${headers} ${sources} )
target_link_libraries (VisionMPC biomimetics lcm qpOASES)

