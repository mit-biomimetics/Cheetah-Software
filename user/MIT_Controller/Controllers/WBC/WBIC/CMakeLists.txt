FILE(GLOB_RECURSE headers *.hpp)
FILE(GLOB_RECURSE sources *.cpp)


add_library(WBIC SHARED ${sources} ${headers})
target_link_libraries(WBIC Goldfarb_Optimizer)
#target_link_libraries(WBIC JCQP)
