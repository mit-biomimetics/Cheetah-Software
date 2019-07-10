FILE(GLOB_RECURSE headers WB*.hpp)
FILE(GLOB_RECURSE sources WB*.cpp)


add_library(WBIC SHARED ${sources} ${headers})
target_link_libraries(WBIC Goldfarb_Optimizer)
target_link_libraries(WBIC JCQP)
