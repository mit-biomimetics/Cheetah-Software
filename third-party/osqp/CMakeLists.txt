# Minimum version required
cmake_minimum_required (VERSION 3.2)

# Project name
project (osqp)

# Export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)

# Some non-standard CMake modules
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/configure/cmake)
INCLUDE(FindPythonModule)
INCLUDE(Utils)
# include(FindMKL)  # Find MKL module


# Detect operating system
# ----------------------------------------------
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(IS_LINUX ON)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(IS_MAC ON)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(IS_WINDOWS ON)
endif()



# Set options
# ----------------------------------------------

# Are unittests generated?
option (UNITTESTS "Enable unittests generation" OFF)

# Is the code generated for embedded platforms?
#   1 :   Yes. Matrix update not allowed.
#   2 :   Yes. Matrix update allowed.
if (NOT DEFINED EMBEDDED)
    message(STATUS "Embedded is OFF")
else()
    message(STATUS "Embedded is ${EMBEDDED}")
endif()

# Is printing enabled?
option (PRINTING "Enable solver printing" ON)
if (DEFINED EMBEDDED)
    message(STATUS "Disabling printing for embedded")
    set(PRINTING OFF)
endif()
message(STATUS "Printing is ${PRINTING}")


# Is profiling enabled?
option (PROFILING "Enable solver profiling (timing)" ON)
if (DEFINED EMBEDDED)
    message(STATUS "Disabling profiling for embedded")
    set(PROFILING OFF)
endif()
message(STATUS "Profiling is ${PROFILING}")

# Is user interrupt enabled?
option (CTRLC "Enable user interrupt (Ctrl-C)" ON)
if (DEFINED EMBEDDED)
    message(STATUS "Disabling user interrupt for embedded")
    set(CTRLC OFF)
endif()
message(STATUS "User interrupt is ${CTRLC}")

# Use floats instead of integers
option (DFLOAT "Use float numbers instead of doubles" OFF)
message(STATUS "Floats are ${DFLOAT}")

# Use long integers for indexing
option (DLONG "Use long integers (64bit) for indexing" ON)
if (NOT (CMAKE_SIZEOF_VOID_P EQUAL 8))
	message(STATUS "Disabling long integers (64bit) on 32bit machine")
	set(DLONG OFF)
endif()
message(STATUS "Long integers (64bit) are ${DLONG}")

# Add code coverage
option (COVERAGE "Perform code coverage" OFF)
message(STATUS "Code coverage is ${COVERAGE}")


# Linear solvers dependencies
# ---------------------------------------------
option (ENABLE_MKL_PARDISO "Enable MKL Pardiso solver" ON)
if (DFLOAT)
	message(STATUS "Disabling MKL Pardiso Solver with floats")
	set(ENABLE_MKL_PARDISO OFF)
endif()
message(STATUS "MKL Pardiso: ${ENABLE_MKL_PARDISO}")


# Generate header file with the global options
# ---------------------------------------------
# modified 2/2/19 by Jared: this doesn't seem to work correctly when osqp is in a subdirectory
# to modify osqp settings, you must modify third-party/osqp/include/osqp_configure.h yourself!
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure/osqp_configure.h.in
#               ${CMAKE_CURRENT_SOURCE_DIR}/include/osqp_configure.h
#               NEWLINE_STYLE LF)

# Set Compiler flags
# ----------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC

# Add specific flag for debug
set(CMAKE_C_FLAGS_DEBUG "-DDDEBUG")

if (NOT MSVC)

    if (COVERAGE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
	if(FORTRAN)
		set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} --coverage")
	endif(FORTRAN)
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
    set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lm")      # Include math
    # Include real time library in linux
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lrt -ldl")
    endif()
endif (NOT MSVC)

# Set sources
# ----------------------------------------------
if(NOT DEFINED EMBEDDED)
    set(
        osqp_src
        src/auxil.c
        src/cs.c
        src/ctrlc.c
        src/kkt.c
        src/lin_alg.c
        src/lin_sys.c
        src/osqp.c
        src/polish.c
        src/proj.c
        src/scaling.c
        src/util.c
    )

    set(
        osqp_headers
        include/auxil.h
        include/constants.h
        include/cs.h
        include/ctrlc.h
        include/glob_opts.h
        include/kkt.h
        include/lin_alg.h
        include/lin_sys.h
        include/osqp.h
        include/osqp_configure.h
        include/types.h
        include/polish.h
        include/proj.h
        include/scaling.h
        include/util.h
    )
else()
    # If embedded is enabled, do not use cs, polish
    set(
        osqp_src
        src/auxil.c
        src/lin_alg.c
        src/lin_sys.c
        src/osqp.c
        src/kkt.c
        src/proj.c
        src/scaling.c
        src/util.c
    )

    set(
        osqp_headers
        include/auxil.h
        include/constants.h
        include/glob_opts.h
        include/lin_alg.h
        include/lin_sys.h
        include/kkt.h
        include/osqp.h
        include/osqp_configure.h
        include/types.h
        include/proj.h
        include/scaling.h
        include/util.h
    )
endif()



# if we are building the Python interface, let's look for Python
# and set some options
# -----------------------------------------------------------------
if (PYTHON)

    # Python include directories need to be passed by the python compilation process
    if (NOT PYTHON_INCLUDE_DIRS)
            message( FATAL_ERROR "You need Python include directories to build the Python interface" )
    endif (NOT PYTHON_INCLUDE_DIRS)

    # Include directories for Python headers
    include_directories(${PYTHON_INCLUDE_DIRS})

    # Pass PYTHON flag to C compiler
    add_definitions(-DPYTHON)

    if (UNITTESTS)
        # Disable unittests
        message(STATUS "Disabling UNITTESTS because we are building Python interface")
        set(UNITTESTS OFF)
    endif (UNITTESTS)

endif (PYTHON)


# if we are building the Matlab interface, let's look for Matlab
# and set some options
# -----------------------------------------------------------------
if (MATLAB)

    find_package(Matlab)

    if (NOT Matlab_FOUND)
        message( FATAL_ERROR "You need Matlab libraries to build the Matlab interface" )
    endif (NOT Matlab_FOUND)

    # Include directories for Matlab headers
    include_directories(${Matlab_INCLUDE_DIRS})

    message(STATUS "Matlab root is " ${Matlab_ROOT_DIR})

    # Pass MATLAB flag to C compiler
    add_definitions(-DMATLAB)

    # Insist on the pre 2018 complex data API
    # so that mxGetPr will work correctly
    add_definitions(-DMATLAB_MEXSRC_RELEASE=R2017b)

    message(STATUS "Using Matlab pre-2018a API for mxGetPr compatibility")

    if (UNITTESTS)
        # Disable unittests
        message(STATUS "Disabling UNITTESTS because we are building Matlab interface")
        set(UNITTESTS OFF)
    endif (UNITTESTS)

endif (MATLAB)

# if we are building the R interface, let's look for R
# and set some options
# -----------------------------------------------------------------
if (R_LANG)

    message(STATUS "We are building the R interface")

    # Look for R libraries
    find_package(R)

    if (NOT R_FOUND)
            message( FATAL_ERROR "You need R libraries to build the R interface" )
    endif (NOT R_FOUND)

    message(STATUS "R exec is: " ${R_EXEC})
    message(STATUS "R root dir is: " ${R_ROOT_DIR})
    message(STATUS "R includes are in: " ${R_INCLUDE_DIRS})

    # Include directories for R headers
    include_directories(${R_INCLUDE_DIRS})

    # Pass R_LANG flag to C compiler
    add_definitions(-DR_LANG)

    if (UNITTESTS)
        # Disable unittests
        message(STATUS "Disabling UNITTESTS because we are building the R interface")
        set(UNITTESTS OFF)
    endif (UNITTESTS)

endif (R_LANG)


# Create Static Library
# ----------------------------------------------

# Add linear system solvers cumulative library
add_subdirectory(lin_sys)

# Static library
#add_library (osqpstatic STATIC ${osqp_src} ${osqp_headers} ${linsys_solvers})
add_library (osqpstatic SHARED ${osqp_src} ${osqp_headers} ${linsys_solvers})
# Give same name to static library output
set_target_properties(osqpstatic PROPERTIES OUTPUT_NAME osqp)

# Include directories for linear system solvers
target_include_directories(osqpstatic PRIVATE ${linsys_solvers_includes})

# Declare include directories for the cmake exported target
target_include_directories(osqpstatic
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/osqp>")

# Install Static Library
# ----------------------------------------------

include(GNUInstallDirs)

install(TARGETS osqpstatic
        EXPORT  ${PROJECT_NAME}
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")


# Install Headers
# ----------------------------------------------

install(FILES ${osqp_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/osqp")



if (MATLAB)
target_link_libraries (osqpstatic ${Matlab_LIBRARIES})
endif (MATLAB)

# If we are building Python/Matlab/R interface:
#   - do not build shared library
#   - do not build demo
if (NOT PYTHON AND NOT MATLAB AND NOT R_LANG)
    # Create osqp shared library
    # NB: Add all the linear system solvers here
    add_library (osqp SHARED ${osqp_src} ${osqp_headers} ${linsys_solvers})

    # Include directories for linear system solvers
    target_include_directories(osqp PRIVATE ${linsys_solvers_includes})

    # Declare include directories for the cmake exported target
    target_include_directories(osqp
                               PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                      "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/osqp>")

    # Install osqp shared library
    install(TARGETS osqp
            EXPORT  ${PROJECT_NAME}
            LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")

    # Create demo executable (linked to static library)
    add_executable (osqp_demo ${PROJECT_SOURCE_DIR}/examples/osqp_demo.c)
    target_link_libraries (osqp_demo osqpstatic)

endif (NOT PYTHON AND NOT MATLAB AND NOT R_LANG)

# Create CMake packages for the build directory
# ----------------------------------------------
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/osqp-targets.cmake"
  NAMESPACE osqp::)

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/osqp-targets.cmake\")\n")
endif()


# Create CMake packages for the install directory
# ----------------------------------------------

set(ConfigPackageLocation lib/cmake/osqp)

install(EXPORT ${PROJECT_NAME}
        FILE osqp-targets.cmake
        NAMESPACE osqp::
        DESTINATION ${ConfigPackageLocation})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake
        DESTINATION ${ConfigPackageLocation})



# Add uninstall command
# ----------------------------------------------
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/configure/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()



# Add testing
# ----------------------------------------------
# Add custom command to generate tests
if (UNITTESTS)
    find_package(PythonInterp)
    if(NOT PYTHONINTERP_FOUND)
        message( FATAL_ERROR "You need python installed to generate unittests. If you do not want to compile the unittests pass -DUNITTESTS=OFF to cmake." )
    endif()

    INCLUDE(FindPythonModule)
    find_python_module(numpy)
    IF(NOT NUMPY_FOUND)
        message( FATAL_ERROR "You need numpy python module installed to generate unittests. If you do not want to compile the unittests pass -DUNITTESTS=OFF to cmake." )
    ENDIF()

    find_python_module(scipy)
    # Check scipy version for sparse.random functionalities
    IF((NOT SCIPY_FOUND) OR (SCIPY_VERSION VERSION_LESS 0.17.0))
        message( FATAL_ERROR "You need scipy python module installed to generate unittests. If you do not want to compile the unittests pass -DUNITTESTS=OFF to cmake." )
    ENDIF()

    find_python_module(__future__)
    IF(NOT __FUTURE___FOUND)
        message( FATAL_ERROR "You need future python module installed to generate unittests. If you do not want to compile the unittests pass -DUNITTESTS=OFF to cmake." )
    ENDIF()

    # Add test_headers and codegen_test_headers
    add_subdirectory(tests)

    # Generating tests.stamp so that the test data are not always generated
    # set(data_timestamp ${PROJECT_SOURCE_DIR}/tests/tests_data.stamp)
    add_custom_command(
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
        COMMAND ${PYTHON_EXECUTABLE} generate_tests_data.py
        DEPENDS ${PROJECT_SOURCE_DIR}/tests/generate_tests_data.py
        OUTPUT ${codegen_test_headers}
        COMMENT "Generating unittests data files using Python"
    )

    # Direct linear solver testing
    include_directories(tests)
    add_executable(osqp_tester
                ${PROJECT_SOURCE_DIR}/tests/osqp_tester.c
                ${PROJECT_SOURCE_DIR}/tests/osqp_tester.h
		${PROJECT_SOURCE_DIR}/tests/minunit.h
                ${test_headers}
                ${codegen_test_headers})
    target_link_libraries (osqp_tester osqpstatic)

    # Add testing
    include(CTest)
    enable_testing()
    add_test(NAME tester COMMAND $<TARGET_FILE:osqp_tester>)
endif()
