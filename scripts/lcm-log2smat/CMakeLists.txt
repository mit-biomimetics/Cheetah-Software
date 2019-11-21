cmake_minimum_required(VERSION 2.6.0)

# pull in the pods macros. See cmake/pods.cmake for documentation
set(POD_NAME lcm-log2smat)
include(cmake/pods.cmake)

# require python
find_package(PythonInterp REQUIRED)

# install all python files in the python/ subdirectory
pods_install_python_packages(${CMAKE_CURRENT_SOURCE_DIR}/python)

# install a script "lcm-log2smat" that runs the lcmlog2smat.log_to_smat python
# module. This script gets installed to ${CMAKE_INSTALL_PREFIX}/bin/lcm-log2smat
# and automatically sets the correct python path.
pods_install_python_script(lcm-log2smat lcmlog2smat.log_to_smat)
