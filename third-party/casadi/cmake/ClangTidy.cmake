function(add_clang_tidy_target TARGET_NAME SOURCES_LIST0)

  list(REMOVE_DUPLICATES SOURCES_LIST0)
  list(SORT SOURCES_LIST0)

  # filtering out unwanted files
  set(SOURCES_LIST)
  foreach(item ${SOURCES_LIST0})
    string(REGEX MATCH ".*\\.hpp" dummy ${item})
    if(NOT dummy)
      string(REGEX MATCH ".*meta\\.cpp" dummy ${item})
      if(NOT dummy)
        list(APPEND SOURCES_LIST ${item})
      endif()
    endif()
  endforeach()

  add_custom_target(clang_tidy_${TARGET_NAME}
    COMMAND clang-tidy
            -p "${PROJECT_BINARY_DIR}/compile_commands.json"
            -warnings-as-errors=*
            "-header-filter=${PROJECT_SOURCE_DIR}/casadi/.*"
            ${SOURCES_LIST}
    DEPENDS ${SOURCES_LIST}
    COMMENT "Clang_tidy ${TARGET_NAME}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  add_custom_target(clang_tidy_fix_${TARGET_NAME}
    COMMAND clang-tidy
            -fix
            -p "${PROJECT_BINARY_DIR}/compile_commands.json"
            "-header-filter=${PROJECT_SOURCE_DIR}/casadi/.*"
            ${SOURCES_LIST}
    DEPENDS ${SOURCES_LIST}
    COMMENT "Clang_tidy ${TARGET_NAME}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)
endfunction()
