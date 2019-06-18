find_package(PythonInterp)

function(add_spell_check_target TARGET_NAME SOURCES_LIST)# PROJECT)

  if(NOT PYTHONINTERP_FOUND)
    return()
  endif()


  list(REMOVE_DUPLICATES SOURCES_LIST)
  list(SORT SOURCES_LIST)

  add_custom_target(spell_${TARGET_NAME}
    COMMAND "${CMAKE_COMMAND}" -E chdir
            "${CMAKE_SOURCE_DIR}/misc"
            "${PYTHON_EXECUTABLE}"
            "spellcheck.py"
            ${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES_LIST}
    DEPENDS ${SOURCES_LIST}
    COMMENT "Spellchecking ${TARGET_NAME}"
    VERBATIM)

endfunction()
