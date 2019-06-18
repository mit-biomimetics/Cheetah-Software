cmake_minimum_required(VERSION 2.8.6)
project(casadi C CXX)

set(CMAKE_CXX_FLAGS "-O3 -Wno-error")
include_directories("${PROJECT_BINARY_DIR}")

set(WITH_IPOPT ON)
set(CASADI_MAJOR_VERSION 3)
set(CASADI_MINOR_VERSION 4)
set(CASADI_PATCH_VERSION 5)
set(CASADI_IS_RELEASE 1)

set(CASADI_VERSION ${CASADI_MAJOR_VERSION}.${CASADI_MINOR_VERSION}.${CASADI_PATCH_VERSION})
if(NOT CASADI_IS_RELEASE)
  set(CASADI_VERSION "${CASADI_VERSION}+")
endif()
set(CASADI_VERSION "${CASADI_VERSION}")
set(PACKAGE_NAME "CasADi")
set(PACKAGE_VERSION ${CASADI_VERSION})
set(PACKAGE_SO_VERSION "3.5")
set(PACKAGE_DESCRIPTION "CasADi is a symbolic framework for automatic differentation and numeric optimization")
set(PACKAGE_AUTHOR "Joel Andersson and Joris Gillis")
set(PACKAGE_MAINTAINER "Joel Andersson and Joris Gillis")
set(PACKAGE_URL "http://casadi.org")

#############################################################
########################## Version ##########################
#############################################################
# For systems with git installed, find out revsion and description.
#execute_process(COMMAND git rev-parse HEAD
  #OUTPUT_VARIABLE git_revision
  #OUTPUT_STRIP_TRAILING_WHITESPACE
  #ERROR_QUIET
#)
#execute_process(COMMAND git describe --first-parent HEAD
  #OUTPUT_VARIABLE git_describe
  #OUTPUT_STRIP_TRAILING_WHITESPACE
  #ERROR_QUIET
#)

#string(REPLACE "+-" "+" git_describe "${git_describe}")
#string(REPLACE "-g" "." git_describe "${git_describe}")
#set(PACKAGE_VERSION_FULL "${git_describe}")

# overwrite version numbers for official release
#if((PACKAGE_VERSION_FULL STREQUAL "") OR (NOT PACKAGE_VERSION MATCHES "\\+"))
  #set(PACKAGE_VERSION_FULL "${PACKAGE_VERSION}")
#endif()

# We use the FeatureSummary module to get a more readable CMake output
include(FeatureSummary)
if(NOT COMMAND add_feature_info)
  function(add_feature_info _name _enabled _description)
    message(STATUS "Detected: ${_name} ${_description}")
  endfunction()
  function(feature_summary)
  endfunction()
endif()

# Default build type is Release
#### WARNING
#if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  #set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
#endif()

########################################################################
############################# Fortran ##################################
########################################################################
option(FORTRAN_REQUIRED "Make fortran a required dependency" ON)

# http://public.kitware.com/Bug/print_bug_page.php?bug_id=9220
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(language_support)
workaround_9220(Fortran Fortran_language_works)
if(FORTRAN_REQUIRED)
  enable_language(Fortran)
  set(Fortran_language_works ON)
else()
  if(Fortran_language_works)
    enable_language(Fortran OPTIONAL)
  endif()
endif()

########################################################################
############################# C++11 ####################################
########################################################################
# Check if the compiler supports C++11
foreach(CXX11FLAG_TEST "" "-std=c++11" "-std=c++0x")
  if(NOT USE_CXX11)
    try_compile(USE_CXX11 # Try to compile test with the proposed flag
      ${CMAKE_BINARY_DIR}
      ${CMAKE_MODULE_PATH}/check_cxx11.cpp
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG_TEST}
      OUTPUT_VARIABLE OUTPUT)

    if(USE_CXX11) # Save flag, if successful
      set(CXX11FLAG ${CXX11FLAG_TEST} CACHE INTERNAL "C++11 compile flag")
    endif()
  endif()
endforeach()

# If found, set flag
if(USE_CXX11)
  if(CXX11FLAG)
    message(STATUS "Flag needed for enabling C++11 features: ${CXX11FLAG}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11FLAG}")
  else()
    message(STATUS "No flag needed for enabling C++11 features.")
  endif()
  add_definitions(-DUSE_CXX11)
else()
  message(FATAL_ERROR "Cannot enable C++11 features.")
endif()

if(CXX11FLAG)
  try_compile(HAS_COPYSIGN
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_copysign.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_COPYSIGN)
    add_definitions(-DHAS_COPYSIGN)
  endif()

  try_compile(HAS_ERF
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erf.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERF)
    add_definitions(-DHAS_ERF)
  endif()

  try_compile(HAS_ERFINV
    ${CMAKE_BINARY_DIR}
    ${CMAKE_MODULE_PATH}/check_erfinv.cpp
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CXX11FLAG}
    OUTPUT_VARIABLE OUTPUT)
  if(HAS_ERFINV)
    add_definitions(-DHAS_ERFINV)
  endif()
endif()

# Check if mkstemps is available
if(WIN32)
  add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
else()
  file(WRITE ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
      "#include <cstdlib>\n"
      "#include <unistd.h>\n"
      "int main() {\n"
      "  char tmp[] = \"mkstempstestXXXXXX.so\";\n"
      "  int flag = mkstemps(tmp, 3);\n"
      "  return flag;\n"
      "}\n")
  try_compile(HAVE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_MKSTEMPS_OUTPUT)
  file(WRITE ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
      "#include <sys/stat.h>\n"
      "#include <fcntl.h>\n"
      "#include <errno.h>\n"
      "int main() {\n"
      "  int fd = open(\"test.c\", O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);\n"
      "}\n")
  try_compile(HAVE_SIMPLE_MKSTEMPS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/simple_mkstemps_test.cxx
              OUTPUT_VARIABLE HAVE_SIMPLE_MKSTEMPS_OUTPUT)
  if(HAVE_MKSTEMPS)
    add_definitions(-DHAVE_MKSTEMPS)
  else()
    if(HAVE_SIMPLE_MKSTEMPS)
      add_definitions(-DHAVE_SIMPLE_MKSTEMPS)
    else()
      message(WARNING "mkstemp is not available, falling back on tmpnam. "
                      "Test compilation returned the following:\n"
                      "${HAVE_MKSTEMPS_OUTPUT}"
                      "Test compilation returned the following:\n"
                      "${HAVE_SIMPLE_MKSTEMPS_OUTPUT}")
   endif()
  endif()
endif()

############################################################
####################### Policies ###########################
############################################################
# Enforce cmake2.4+ style behaviour for
# mixed filename/full path linker directives
cmake_policy(SET CMP0003 NEW)

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

# cmake_policy avoid escaping in add_definitions
cmake_policy(SET CMP0005 OLD)

# cmake_policy use MACOSX_RPATH=ON by default
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

# cmake_policy how to interpret string arguments
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

############################################################
################# WITH_COMMON option #######################
############################################################
option(WITH_COMMON "Build common, stable interfaces, if found" ON)
if(WITH_COMMON)
  # Look for LAPACK
  find_package(LAPACK QUIET)
  set(WITH_LAPACK_DEF ${LAPACK_FOUND})
  # qpOASES and blockSQP built if LAPACK available
  set(WITH_QPOASES_DEF ${LAPACK_FOUND})
  set(WITH_BLOCKSQP_DEF ${LAPACK_FOUND})
  # Look for IPOPT
  find_package(IPOPT QUIET)
  set(WITH_IPOPT_DEF ${IPOPT_FOUND})
else()
  # Disabled by default
  set(WITH_LAPACK_DEF OFF)
  set(WITH_QPOASES_DEF OFF)
  set(WITH_BLOCKSQP_DEF OFF)
  set(WITH_IPOPT_DEF OFF)
endif()

############################################################
################# Misc options #############################
############################################################
option(ENABLE_EXPORT_ALL "Export all symbols to a shared library" OFF)
option(ENABLE_STATIC "Compile static libraries (.a on Linux)" OFF)
option(ENABLE_SHARED "Compile shared libraries (.so on Linux)" ON)
option(SWIG_EXPORT "Export SWIG" OFF)
option(SWIG_IMPORT "Export SWIG" OFF)

# enabling OpenMP support if requested
option(WITH_OPENMP "Compile with parallelization support" OFF)
if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DWITH_OPENMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DWITH_OPENMP")
  endif()
endif()

# enabling Threads support if requested
option(WITH_THREAD "Compile with parallelization support" ON)
option(WITH_THREAD_MINGW "For MINGW, use mingw-std-threads" OFF)
if(WITH_THREAD)
  add_definitions(-DCASADI_WITH_THREAD)
endif()
if(MINGW AND WITH_THREAD_MINGW)
  add_definitions(-DCASADI_WITH_THREAD_MINGW)
else()
  if(WITH_THREAD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()


# OpenCL
option(WITH_OPENCL "Compile with OpenCL support (experimental)" OFF)
if(WITH_OPENCL)
  # Core depends on OpenCL for GPU calculations
  find_package(OpenCL REQUIRED)
  add_definitions(-DWITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()
add_feature_info(opencl-support WITH_OPENCL "Enable just-in-time compiliation to CPUs and GPUs with OpenCL.")

# Enable: RTLD_DEEPBIND
option(WITH_DEEPBIND "Load plugins with RTLD_DEEPBIND (can be used to resolve conflicting libraries in e.g. MATLAB)" ON)
if(WITH_DEEPBIND)
  add_definitions(-DWITH_DEEPBIND)
endif()

# Suppress warnings for refcounting failures (at program closure)
option(WITH_REFCOUNT_WARNINGS "Issue warnings upon reference counting failure" OFF)
if(WITH_REFCOUNT_WARNINGS)
  add_definitions(-DWITH_REFCOUNT_WARNINGS)
endif()

# Have an so version?
option(WITH_SO_VERSION "Use an so version for the library (version suffix) when applicable" ON)

############################################################
################# Installation location ####################
############################################################
# Flat hierarchy if self-contained and on Windows
option(WITH_SELFCONTAINED "Make the install directory self-contained" OFF)
if(WITH_SELFCONTAINED OR (WIN32 AND NOT CYGWIN))
  set(DEF_LIB_PREFIX casadi)
  set(DEF_BIN_PREFIX casadi)
  set(DEF_CMAKE_PREFIX casadi/cmake)
  set(DEF_PYTHON_PREFIX) # Note: Files installed in a casadi subdirectory
  set(DEF_MATLAB_PREFIX casadi)
  set(DEF_INCLUDE_PREFIX casadi/include)
else()
  set(DEF_LIB_PREFIX lib)
  set(DEF_BIN_PREFIX bin)
  set(DEF_CMAKE_PREFIX lib/cmake/casadi)
  set(DEF_PYTHON_PREFIX python)
  if (WITH_PYTHON)
    set(PYTHON_VERSION_REQ "")
    if (WITH_PYTHON3)
      set(PYTHON_VERSION_REQ "3")
    endif()
    find_package(PythonInterp ${PYTHON_VERSION_REQ} REQUIRED)
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; from distutils.sysconfig import get_python_lib; print(get_python_lib())"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(DEF_PYTHON_PREFIX ${PYTHON_SITE_PACKAGES})
  endif()
  set(DEF_MATLAB_PREFIX matlab)
  set(DEF_INCLUDE_PREFIX include)
endif()

# Installation directories are cache entries that can be modified
set(LIB_PREFIX ${DEF_LIB_PREFIX} CACHE PATH "Installation directory for libraries")
set(BIN_PREFIX ${DEF_BIN_PREFIX} CACHE PATH "Installation directory for executables")
set(INCLUDE_PREFIX ${DEF_INCLUDE_PREFIX} CACHE PATH "Installation directory for header files")
set(CMAKE_PREFIX ${DEF_CMAKE_PREFIX} CACHE PATH "Installation directory for CMake files")
set(PYTHON_PREFIX ${DEF_PYTHON_PREFIX} CACHE PATH "Installation directory for Python files")
set(MATLAB_PREFIX ${DEF_MATLAB_PREFIX} CACHE PATH "Installation directory for MATLAB/Octave files")
add_definitions(-DPLUGIN_EXTRA_SEARCH_PATH="\\"${CMAKE_BINARY_DIR}/third-party/casadi/lib\\"")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE PYTHON MATLAB)
  set(var ${p}_PREFIX)
  if(${var})
    if(NOT IS_ABSOLUTE "${${var}}")
      # Turn relative path into absolute
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
  else()
    # Flat installation hierarchy
    set(${var} "${CMAKE_INSTALL_PREFIX}")
  endif()
endforeach()

# Set rpath
if(WITH_SELFCONTAINED)
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
  endif()
else()
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIB_PREFIX}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${LIB_PREFIX}")
  endif()
endif()

################################################################################
############################### compiler warnings ##############################
################################################################################
# add compiler warnings
option(WITH_EXTRA_WARNINGS "Add extra warnings (like -Wall -Wextra)" OFF)
option(WITH_WERROR "Add -Werror" OFF)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-mismatched-tags -Wno-tautological-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  if(WITH_EXTRA_WARNINGS)
    set(EXTRA_CASADI_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare")
    set(EXTRA_CASADI_C_FLAGS "-Wall -Wextra -Wformat-security -Wno-unused-parameter -Wno-maybe-uninitialized -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value")
  endif()
  if(WITH_WERROR)
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Werror")
    set(EXTRA_CASADI_C_FLAGS "${EXTRA_CASADI_C_FLAGS} -Werror")
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # I don't know the MSVC warnings, so do nothing for now
endif()

#######################################################################
############################### Extra runtime checks ##################
#######################################################################
option(WITH_EXTRA_CHECKS "Compile in extra runtime checks (for developers)" OFF)
if(WITH_EXTRA_CHECKS)
   add_definitions(-DWITH_EXTRA_CHECKS)
endif()


#######################################################################
############################# -fPIC ###################################
#######################################################################
# the following commands are needed to fix a problem with the libraries
# for linux 64 bits
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  message(STATUS "x86_64 architecture detected - setting flag -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
endif()


include_directories(.)
include_directories(${CMAKE_BINARY_DIR})

#######################################################################
########################### code coverage #############################
#######################################################################
option(WITH_COVERAGE "Create coverage report" OFF)
if(WITH_COVERAGE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS  " -lgcov -fprofile-arcs --coverage ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(MINGW)
  # Circumventing a bug in MinGW g++ v4.7.2, evoked by 752fa89355ffa
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-ipa-cp-clone")
endif()

# Tentative fix for https://github.com/casadi/casadi/issues/1665
if(MINGW OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fkeep-inline-functions")
endif()

# Plugin version
add_definitions(-DCASADI_VERSION=31)

# Necessary for Visual C++
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

# Dynamic loading (needed for loading plugins and external functions)
option(WITH_DL "Enable dynamic loading of functions" ON)
if(WITH_DL)
  add_definitions(-DWITH_DL)
endif()
add_feature_info(dynamic-loading WITH_DL "Compile with support for dynamic loading of generated functions (needed for ExternalFunction)")

# Include support for deprecated features (to be removed in the next release)
option(WITH_DEPRECATED_FEATURES "Compile with syntax that is scheduled to be deprecated" ON)
if (WITH_DEPRECATED_FEATURES)
  add_definitions(-DWITH_DEPRECATED_FEATURES)
endif()

include_directories(.)
include_directories(${CMAKE_BINARY_DIR})


#######################################################################
################# 3rd party interface #################################
#######################################################################
# Sundials: Solve initial-value problems in ODE and DAE
option(WITH_BUILD_SUNDIALS "Compile the included source code for Sundials 2.6.1" ON)
option(WITH_SUNDIALS "Compile the interface to Sundials" ON)
if(WITH_SUNDIALS)
  if(WITH_BUILD_SUNDIALS)
    # build the included sundials
    set(SUNDIALS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/casadi-sundials)
    set(SUNDIALS_LIBRARIES casadi_sundials)
  else()
    # try to find system sundials
    add_definitions(-DWITH_SYSTEM_SUNDIALS)
    find_package(SUNDIALS REQUIRED)
  endif()
endif()
add_feature_info(sundials-interface WITH_SUNDIALS "Interface to the ODE/DAE integrator suite SUNDIALS.")

# CSparse: A sparse direct linear solver (included in the source)
option(WITH_BUILD_CSPARSE "Compile the included source code for CSparse" ON)
option(WITH_CSPARSE "Compile the interface to CSparse" ON)
if(WITH_CSPARSE)
  if(WITH_BUILD_CSPARSE)
    # build the included csparse
    set(CSPARSE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/CSparse/Include)
    set(CSPARSE_LIBRARIES casadi_csparse)
  else()
    # try to find system csparse
    find_package(CSPARSE REQUIRED)
  endif()
endif()
add_feature_info(csparse-interface WITH_CSPARSE "Interface to the sparse direct linear solver CSparse.")

# BLASFEO
option(WITH_BLASFEO "Compile the interface to BLASFEO" OFF)
option(WITH_BUILD_BLASFEO "Compile the the included source code for BLASFEO" ON)
if(WITH_BLASFEO)
  if(WITH_BUILD_BLASFEO)
    if(EXISTS ${PROJECT_SOURCE_DIR}/external_packages/blasfeo/src/include)
      set(BLASFEO_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/blasfeo/src/include)
      set(BLASFEO_LIBRARIES)
    else()
      message(FATAL_ERROR "Blasfeo directory not found. If you want blasfeo, do git submodule init; git submodule update.")
    endif()
  else()
    # try to find system blasfeo
    find_package(BLASFEO REQUIRED)
    set(BLASFEO_DLOPEN OFF)
  endif()
endif()
add_feature_info(blasfeo-interface WITH_BLASFEO "Interface to blasfeo.")

# HPMPC
option(WITH_HPMPC "Compile the interface to HPMPC" OFF)
option(WITH_BUILD_HPMPC "Compile the the included source code for HPMPC" ON)
if(WITH_HPMPC)
  if(NOT WITH_BLASFEO)
    message(FATAL_ERROR "WITH_HPMPC requires WITH_BLASFEO")
  endif()
  if(WITH_BUILD_HPMPC)
    if(EXISTS ${PROJECT_SOURCE_DIR}/external_packages/hpmpc/src/include)
      set(HPMPC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/hpmpc/src/include)
      set(HPMPC_LIBRARIES)
      set(HPMPC_DLOPEN ON)
    else()
      message(FATAL_ERROR "Hpmpc directory not found. If you want hpmpc, do git submodule init; git submodule update.")
    endif()
  else()
    # try to find system hpmpc
    find_package(HPMPC REQUIRED)
    set(HPMPC_DLOPEN OFF)
  endif()
endif()
add_feature_info(hpmpc-interface WITH_HPMPC "Interface to hpmpc.")

# TinyXML: An XML reader, needed e.g. by the Modelica interface
option(WITH_BUILD_TINYXML "Compile the included TinyXML source code" ON)
option(WITH_TINYXML "Compile the interface to TinyXML" ON)
if(WITH_TINYXML)
  if(WITH_BUILD_TINYXML)
    # build the included tinyxml
    set(TINYXML_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/tinyxml-2.6.2)
    set(TINYXML_LIBRARIES casadi_tinyxml)
  else()
    # try to find system tinyxml
    find_package(TINYXML REQUIRED)
  endif()
endif()
add_feature_info(tinyxml-interface WITH_TINYXML "Interface to the XML parser TinyXML.")

# DSDP
option(WITH_BUILD_DSDP "Compile the the included source code for DSDP" OFF)
option(WITH_DSDP "Compile the interface to DSDP" OFF)
if(WITH_DSDP)
  if(WITH_BUILD_DSDP)
    # build the included dsdp
    set(DSDP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/DSDP/include)
    set(DSDP_LIBRARIES casadi_dsdp)
  else()
    # try to find system dsdp
    find_package(DSDP REQUIRED)
  endif()
endif()
add_feature_info(dsdp-interface WITH_DSDP "Interface to the interior point SDP solver DSDP.")

# CLANG: A just-in-time compiler
option(OLD_LLVM "Use pre-3.5 LLVM" OFF)
option(WITH_CLANG "Compile the interface to clang JIT" OFF)
if(WITH_CLANG)
  # try to find system CLANG
  find_package(CLANG REQUIRED)
endif()
add_feature_info(clang-interface WITH_CLANG "Interface to the Clang JIT compiler.")

# Lapack: Dense linear solvers
option(WITH_LAPACK "Compile the interface to LAPACK" ${WITH_LAPACK_DEF})
if(WITH_LAPACK)
  find_package(LAPACK REQUIRED)
  if(DEFINED ENV{CASADI_EXTRA_LIBRARIES_LAPACK})
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} $ENV{CASADI_EXTRA_LIBRARIES_LAPACK})
  endif()
endif()
add_feature_info(lapack-interface WITH_LAPACK "Interface to LAPACK.")

# qpOASES: An active-set QP solver
option(WITH_QPOASES "Compile the interface to qpOASES (the source code for qpOASES is included)" ${WITH_QPOASES_DEF})
option(WITH_NO_QPOASES_BANNER "Add -D__NO_COPYRIGHT__ to qpOASES definitions" OFF)
if(WITH_QPOASES)
  if (NOT WITH_LAPACK)
     message(FATAL_ERROR "WITH_QPOASES requires WITH_LAPACK.")
  endif()
  set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/qpOASES/include)
  set(QPOASES_LIBRARIES casadi_qpoases)
endif()
add_feature_info(qpoases-interface WITH_QPOASES "Interface to the active-set QP solver qpOASES.")

# blockSQP: A block-structure-exploiting SQP method
option(WITH_BLOCKSQP "Compile the interface to blockSQP (the source code for blockSQP is included)" ${WITH_BLOCKSQP_DEF})
if(WITH_BLOCKSQP)
  if (NOT WITH_QPOASES)
    message(FATAL_ERROR "WITH_BLOCKSQP requires WITH_QPOASES.")
  endif()
endif()
add_feature_info(blocksqp-interface WITH_BLOCKSQP "Interface to the NLP solver blockSQP.")

# IPOPT: A sparse interior point NLP solver with filter line-search
option(WITH_IPOPT "Compile the interface to IPOPT" ${WITH_IPOPT_DEF})
if(WITH_IPOPT)
  find_package(IPOPT REQUIRED)
endif()
add_feature_info(ipopt-interface WITH_IPOPT "Interface to the NLP solver Ipopt.")

# Bonmin: An MINLP solver
option(WITH_BONMIN "Compile the interface to BONMIN" OFF)
if(WITH_BONMIN)
  find_package(BONMIN REQUIRED)
endif()
add_feature_info(bonmin-interface WITH_BONMIN "Interface to the MINLP framework Bonmin.")

# Knitro: An NLP solver
option(WITH_KNITRO "Compile the interface to KNITRO" OFF)
if(WITH_KNITRO)
  find_package(KNITRO REQUIRED)
endif()
add_feature_info(knitro-interface WITH_KNITRO "Interface to the NLP solver KNITRO.")

# CPLEX: A sparse mixed-integer QP solver
option(WITH_CPLEX "Compile the interface to CPLEX" OFF)
option(WITH_CPLEX_SHARED "Use shared cplex libraries" OFF)
if(WITH_CPLEX)
  find_package(CPLEX REQUIRED)
endif()
add_feature_info(cplex-interface WITH_CPLEX "Interface to the QP solver CPLEX.")

# CLP: A sparse linear solver
option(WITH_CLP "Compile the CLP interface" OFF)
if(WITH_CLP)
  find_package(CLP REQUIRED)
endif()
add_feature_info(clp-interface WITH_CLP "Interface to the LP solver CLP.")

# SNOPT: An SQP method
option(WITH_SNOPT "Compile the SNOPT interface if SNOPT libraries can be found" OFF)
option(WITH_SNOPT_FORCE_DUMMY "Compile the SNOPT interface against a dummy library (for Debian distribution; breaks normal SNOPT build)" OFF)
option(WITH_BUILD_SNOPT_INTERFACE "Use the SNOPT C-interface shipped with casadi" ON)
if(WITH_SNOPT)
  if(NOT WITH_SNOPT_FORCE_DUMMY)
    find_package(SNOPT REQUIRED)
  endif()
  if(WITH_BUILD_SNOPT_INTERFACE)
    set(SNOPT_INTERFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_packages/snopt-interface/include)
    set(SNOPT_INTERFACE_LIBRARIES casadi_snopt_interface)
  else()
    find_package(SNOPT_INTERFACE REQUIRED)
  endif()
endif()
add_feature_info(snopt-interface WITH_SNOPT "Interface to the NLP solver SNOPT.")

# HSL: Sparse direct linear solvers
option(WITH_HSL "Enable HSL interface" OFF)
if(WITH_HSL)
  find_package(HSL REQUIRED)
endif()
add_feature_info(hsl-interface WITH_HSL "Interface to HSL.")

# OOQP: A sparse interior point QP solver
option(WITH_OOQP "Enable OOQP interface" OFF)
if(WITH_OOQP)
  find_package(OOQP REQUIRED)
endif()
add_feature_info(ooqp-interface WITH_OOQP "Interface to the QP solver OOQP (requires BLAS and HSL libraries).")

# GUROBI: A sparse mixed-integer QP solver
option(WITH_GUROBI "Compile the interface to GUROBI" OFF)
if(WITH_GUROBI)
  find_package(GUROBI REQUIRED)
endif()
add_feature_info(gurobi-interface WITH_GUROBI "Interface to the (mixed-integer) QP solver GUROBI")

# SQIC: A QP solver
option(WITH_SQIC "Enable SQIC interface" OFF)
if(WITH_SQIC)
  if (NOT Fortran_language_works)
    message(FATAL_ERROR "SQIC requires Fortran support.")
  endif()
  find_package(SQIC REQUIRED)
endif()
add_feature_info(sqic-interface WITH_SQIC "Interface to the QP solver SQIC.")

# AMPL: Interface to the AMPL solver library and interfaced solvers
option(WITH_AMPL "Compile the interface to AMPL" OFF)
add_feature_info(ampl-interface WITH_AMPL "Interface to the AMPL solver library.")

# Slicot: used for expm, dple
option(WITH_SLICOT "Enable SLICOT interface" OFF)
if(WITH_SLICOT)
  if (NOT Fortran_language_works)
    message(FATAL_ERROR "SLICOT requires Fortran support.")
  endif()
  if (NOT WITH_LAPACK)
     message(FATAL_ERROR "WITH_SLICOT requires WITH_LAPACK.")
  endif()
  find_package(SLICOT REQUIRED)
endif()
add_feature_info(slicot-interface WITH_SLICOT "Interface to the controls library SLICOT.")

# Worhp: A sparse SQP method
option(WITH_WORHP "Compile the WORHP interface" OFF)
if(WITH_WORHP)
  find_package(WORHP REQUIRED)
endif()
add_feature_info(worhp-interface WITH_WORHP "Interface to the NLP solver Worhp (requires gfortran, gomp).")

######################################################
##################### paths ##########################
######################################################
if(WIN32)
  # On Windows, build everything in the same directory
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
else()
  # On non-Windows, use separate directories for build
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
endif()
#####################################################
######################### swig ######################
#####################################################
option(WITH_PYTHON "Compile the Python front-end" OFF)
option(WITH_MATLAB "Compile the MATLAB front-end" OFF)
option(WITH_OCTAVE "Compile the Octave front-end (experimental)" OFF)
option(WITH_JSON "Compile the JSON front-end" OFF)
option(WITH_OCTAVE_IMPORT "[DEPRECATED] Copy import.m file to install directory" OFF)
option(WITH_PYTHON3 "Use Python 3.*" OFF)
option(WITH_EXTENDING_CASADI "Compile a demonstration that shows how a project that depends on CasADi can be implemented." OFF)
if(WITH_PYTHON OR WITH_MATLAB OR WITH_OCTAVE OR WITH_JSON)
  add_subdirectory(swig)
endif()

option(WITH_LINT "Compile with support for linting (syntax checking)" OFF)
if(WITH_LINT)
  include(CppLint REQUIRED)
endif()

option(WITH_SPELL "Compile with support for spell-checking" OFF)
if(WITH_SPELL)
  include(SpellCheck REQUIRED)
endif()

option(WITH_CLANG_TIDY "Compile with support for clang-tidy" OFF)
if(WITH_CLANG_TIDY)
  include(ClangTidy REQUIRED)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if(WITH_EXTRA_WARNINGS)
    # -Wall leaks into clang-tidy (through compile_commands.json)
    # silently overruling the 'Checks' field.
    # These warnings come from external_packages
    set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Wno-overloaded-virtual -Wno-deprecated-register")
  endif()
endif()

include(embed_resources)
include(GenerateExportHeader)


#####################################################
################# include the rest ##################
#####################################################
add_subdirectory(external_packages)
add_subdirectory(experimental EXCLUDE_FROM_ALL)
add_subdirectory(misc)
add_subdirectory(casadi)

# Print a summary
feature_summary(WHAT ALL)
feature_summary(WHAT ENABLED_FEATURES VAR feature_list)

string(REPLACE "\n" "\\n" feature_list "${feature_list}")
get_directory_property( compiler_definitions COMPILE_DEFINITIONS )
set(EXTRA_CXX_FLAGS_FROM_DEFS "")
foreach(DEF ${compiler_definitions})
  set(EXTRA_CXX_FLAGS_FROM_DEFS "${EXTRA_CXX_FLAGS_FROM_DEFS} -D${DEF}")
endforeach()

# what is this?
if(EXTRA_CMAKE)
  include(${EXTRA_CMAKE})
endif()

########################################################
########### Generate config files ######################
########################################################
# Add all targets to the build-tree export set
export(TARGETS casadi
  FILE "${PROJECT_BINARY_DIR}/casadi-targets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE CasADi)

# Create the casadi-config.cmake and casadi-config-version.cmake files
file(RELATIVE_PATH REL_INCLUDE_PREFIX "${CMAKE_PREFIX}" "${INCLUDE_PREFIX}")
file(RELATIVE_PATH REL_LIB_PREFIX "${CMAKE_PREFIX}" "${LIB_PREFIX}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "${LIBRARY_OUTPUT_PATH}")
configure_file(casadi-config.cmake.in "${PROJECT_BINARY_DIR}/casadi-config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${CASADI_CMAKE_DIR}/${REL_INCLUDE_PREFIX}")
set(CONF_LIBRARY_DIRS "\${CASADI_CMAKE_DIR}/${REL_LIB_PREFIX}")
configure_file(casadi-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake" @ONLY)
# ... for both
configure_file(casadi-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake" @ONLY)

# Install the casadi-config.cmake and casadi-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/casadi-config.cmake"
  "${PROJECT_BINARY_DIR}/casadi-config-version.cmake"
  DESTINATION "${CMAKE_PREFIX}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT casadi-targets DESTINATION "${CMAKE_PREFIX}" COMPONENT dev)

set(CPACK_PACKAGE_CONTACT "casadi-users@googlegroups.com")
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION_FULL})
include(CPack)
