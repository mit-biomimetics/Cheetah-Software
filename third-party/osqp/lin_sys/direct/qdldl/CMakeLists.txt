# Add qdldl
add_subdirectory(qdldl_sources)


if(NOT DEFINED EMBEDDED)
set(
    amd_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/include/amd_internal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/include/amd.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/include/SuiteSparse_config.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_aat.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_control.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_defaults.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_info.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_order.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_post_tree.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_postorder.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_preprocess.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/amd_valid.c
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/src/SuiteSparse_config.c
)
endif()


set(qdldl_interface_includes
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/amd/include
    ${CMAKE_CURRENT_SOURCE_DIR}/qdldl_sources/include
)

set(qdldl_interface_src
    ${amd_sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/qdldl_interface.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qdldl_interface.c
)

# Create object library for linear system solver interface
add_library(linsys_qdldl OBJECT ${qdldl_interface_src})
target_include_directories(linsys_qdldl PRIVATE ${qdldl_interface_includes} ${PROJECT_SOURCE_DIR}/include)
