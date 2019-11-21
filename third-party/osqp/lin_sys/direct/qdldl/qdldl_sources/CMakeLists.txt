# Minimum version required
cmake_minimum_required (VERSION 3.2)

# Project name
project (qdldl)

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)

# Some non-standard CMake modules
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/configure/cmake)

# Export compilation commands for IDEs and autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
# ----------------------------------------------
# Use floats instead of doubles
set(DFLOAT OFF CACHE BOOL "Use float numbers instead of doubles")
message(STATUS "Floats are ${DFLOAT}")

# Use long integers for indexing
set(DLONG ON CACHE BOOL "Use long integers (64bit) for indexing")
if (NOT (CMAKE_SIZEOF_VOID_P EQUAL 8))
	message(STATUS "Disabling long integers (64bit) on 32bit machine")
	set(DLONG OFF)
endif()
message(STATUS "Long integers (64bit) are ${DLONG}")


# Set Compiler flags
# ----------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC


# Add compiler options if we are not on windows
if (NOT MSVC)

    if (COVERAGE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
	if(FORTRAN)
		set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} --coverage")
	endif(FORTRAN)
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
endif (NOT MSVC)

# Generate header file with the global options
# ---------------------------------------------

# numeric types
if(DFLOAT)
  set(QDLDL_FLOAT_TYPE "float")
else()
	set(QDLDL_FLOAT_TYPE "double")
endif()

if(DLONG)
  set(QDLDL_INT_TYPE "long long")
else()
	set(QDLDL_INT_TYPE "int")
endif()

#boolean type is always unsigned char
#for now, since _Bool does not exist in
#C89 and we want to avoid interoperability
#problems when calling QDLDL from C++
set(QDLDL_BOOL_TYPE "unsigned char")

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure/qdldl_types.h.in
#               ${CMAKE_CURRENT_BINARY_DIR}/include/qdldl_types.h
#               NEWLINE_STYLE LF)


# Set sources
# ----------------------------------------------
set(
	qdldl_src
	src/qdldl.c
	)

set(
	qdldl_headers
	include/qdldl.h
	include/qdldl_types.h
	)

# Create object library
# ----------------------------------------------
add_library (qdldlobject OBJECT ${qdldl_src} ${qdldl_headers})
target_include_directories(qdldlobject PRIVATE ${PROJECT_SOURCE_DIR}/include)


# Create Static Library
# ----------------------------------------------

# Static library
add_library (qdldlstatic STATIC ${qdldl_src} ${qdldl_headers})
# Give same name to static library output
set_target_properties(qdldlstatic PROPERTIES OUTPUT_NAME qdldl)

# Declare include directories for the cmake exported target
target_include_directories(qdldlstatic
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/qdldl>")

# Install Static Library
# ----------------------------------------------

include(GNUInstallDirs)

install(TARGETS qdldlstatic
        EXPORT  ${PROJECT_NAME}
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")


# Install Headers
# ----------------------------------------------
install(FILES ${qdldl_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/qdldl")


# Install Shared Library
# ----------------------------------------------
# Create qdldl shared library
add_library (qdldl SHARED ${qdldl_src} ${qdldl_headers})

# Declare include directories for the cmake exported target
target_include_directories(qdldl
	PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/qdldl>")

# Install qdldl shared library
install(TARGETS qdldl
	EXPORT  ${PROJECT_NAME}
	LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Create demo executable (linked to static library)
add_executable (qdldl_example ${PROJECT_SOURCE_DIR}/examples/c/example.c)
target_link_libraries (qdldl_example qdldlstatic)


# Create CMake packages for the build directory
# ----------------------------------------------
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/qdldl-targets.cmake"
  NAMESPACE qdldl::)

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/qdldl-targets.cmake\")\n")
endif()


# Create CMake packages for the install directory
# ----------------------------------------------

set(ConfigPackageLocation lib/cmake/qdldl)

install(EXPORT ${PROJECT_NAME}
        FILE qdldl-targets.cmake
        NAMESPACE qdldl::
        DESTINATION ${ConfigPackageLocation})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/qdldl-config.cmake
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

    # Add test_headers and codegen_test_headers
    add_subdirectory(tests)

    # Direct qdldl solver testing
    add_executable(qdldl_tester
                ${PROJECT_SOURCE_DIR}/tests/qdldl_tester.c ${PROJECT_SOURCE_DIR}/tests/minunit.h
                ${test_headers})
    target_link_libraries (qdldl_tester qdldlstatic)

    # Add testing
    include(CTest)
    enable_testing()
    add_test(NAME tester COMMAND $<TARGET_FILE:qdldl_tester>)
endif()
