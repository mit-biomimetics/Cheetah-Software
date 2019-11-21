cmake_minimum_required(VERSION 3.5)
project(sim)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)


execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/find_qt_path.sh OUTPUT_VARIABLE CMAKE_PREFIX_PATH_1)
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH_1})
# message("CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 OLD)
endif()

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)
find_package(Qt5Core CONFIG REQUIRED)
find_package(Qt5Widgets CONFIG REQUIRED)
get_target_property(QtCore_location Qt5::Core LOCATION)
message(STATUS ${QtCore_location})
find_package(Qt5Gui CONFIG REQUIRED)

set(QOBJECT_HEADERS
        include/Graphics3D.h
        include/GameController.h
        include/SimControlPanel.h
        include/JoystickTest.h)
#        include/main_window.h
#        include/gamepad_monitor.h
#        include/spi_debug.h)

QT5_WRAP_CPP(tst_hdr_moc ${QOBJECT_HEADERS})

include_directories("include/")
include_directories("../common/include/")
include_directories("../lcm-types/cpp")
include_directories("/usr/local/include/lcm/")   # lcm includes
include_directories("../third-party/ParamHandler")
include_directories(${CMAKE_BINARY_DIR})
file(GLOB sources "src/*.cpp")

find_package(Qt5Core COMPONENTS QtGamepad REQUIRED)
find_package(Qt5Gamepad REQUIRED)

set(OpenGL_GL_PREFERENCE LEGACY)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

add_executable(sim ${sources} ${tst_hdr_moc})
target_link_libraries(sim biomimetics pthread lcm inih dynacore_param_handler 
    Qt5::Widgets Qt5::Core Qt5::Gui Qt5::Gamepad ${QT_LIBRARIES} ${OPENGL_LIBRARIES})

if (CMAKE_SYSTEM_NAME MATCHES Linux)
  target_link_libraries(sim rt ${GLUT_glut_LIBRARIES})
  include (CheckIncludeFileCXX)
endif (CMAKE_SYSTEM_NAME MATCHES Linux)

if (APPLE)
  target_link_libraries(sim "-framework GLUT")
endif (APPLE)


