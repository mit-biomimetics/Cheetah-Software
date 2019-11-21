FILE(GLOB_RECURSE headers *.hh)
FILE(GLOB_RECURSE sources *.cc)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb ") # no werror here!


add_library (Goldfarb_Optimizer SHARED ${headers} ${sources})
#add_library (Goldfarb_Optimizer STATIC ${headers} ${sources})

#install(TARGETS Goldfarb_Optimizer DESTINATION "${INSTALL_LIB_DIR}")
#install(FILES ${headers} DESTINATION
    #"${INSTALL_INCLUDE_DIR}/Optimizer/Goldfarb_Optimizer")
