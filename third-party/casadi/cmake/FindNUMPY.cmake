if(NOT NUMPY_PATH)
  find_package(PythonLibs QUIET)
  find_package(PythonInterp QUIET)
  if(PYTHON_EXECUTABLE AND NOT NUMPY_INCLUDE_DIR)
    message(STATUS "Python executable is ${PYTHON_EXECUTABLE}")
    exec_program ("${PYTHON_EXECUTABLE}" 
         ARGS "-c \"import numpy; print(numpy.get_include())\""
         OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
         RETURN_VALUE NUMPY_NOT_FOUND)
    message(STATUS "numpy.get_include() is ${NUMPY_INCLUDE_DIR}")
  endif(PYTHON_EXECUTABLE AND NOT NUMPY_INCLUDE_DIR)
  
  find_path(NUMPY_PATH
  ndarrayobject.h
  HINTS $ENV{NUMPY_INCLUDE}
  PATHS ${PYTHON_INCLUDE_DIRS}/numpy ${PYTHON_INCLUDE_PATH}/numpy ${NUMPY_INCLUDE_DIR}/numpy 
  )
endif()

# Joel: What is the purpose of this??
find_path(NUMPY_PATH_DEP1
endian.h
)

string(REGEX REPLACE "include" "lib" NUMPY_LIBS_INTERMEDIATE "${NUMPY_PATH}")
string(REGEX REPLACE "numpy$" "" NUMPY_LIBS_BASE "${NUMPY_LIBS_INTERMEDIATE}")

string(REGEX REPLACE "include" "lib/pyshared" NUMPY_LIBS_INTERMEDIATE_2 "${NUMPY_PATH}")
string(REGEX REPLACE "numpy$" "" NUMPY_LIBS_BASE_2 "${NUMPY_LIBS_INTERMEDIATE_2}")

string(REGEX REPLACE "include" "lib/pymodules" NUMPY_LIBS_INTERMEDIATE_3 "${NUMPY_PATH}")
string(REGEX REPLACE "numpy$" "" NUMPY_LIBS_BASE_3 "${NUMPY_LIBS_INTERMEDIATE_3}")

find_library(NUMPY_LIBS_1
NAMES multiarray multiarray.so
HINTS $ENV{NUMPY_LIBS}
PATHS ${NUMPY_LIBS_BASE}dist-packages/numpy/core
      ${NUMPY_LIBS_BASE_2}numpy/core
      ${NUMPY_LIBS_BASE_3}numpy/core
      ${NUMPY_INCLUDE_DIR}/..
)

set(NUMPY_LIBS "${NUMPY_LIBS_1}")
if(NOT NUMPY_INCLUDED_DIRS)
  set(NUMPY_INCLUDED_DIRS ${NUMPY_PATH})
  if(NUMPY_PATH_DEP1)
    set(NUMPY_INCLUDED_DIRS "${NUMPY_INCLUDED_DIRS};${NUMPY_PATH_DEP1}")
  endif(NUMPY_PATH_DEP1)
endif()

if(NUMPY_PATH OR NUMPY_INCLUDED_DIRS)
set(NUMPY_FOUND TRUE)
message(STATUS "Numpy path found: ${NUMPY_PATH}")
message(STATUS "Python libs: ${PYTHON_LIBRARIES}")
message(STATUS "Numpy includes: ${NUMPY_INCLUDED_DIRS}")
message(STATUS "Numpy libs: ${NUMPY_LIBS}")
else()
message(STATUS "Numpy not found")
endif()
