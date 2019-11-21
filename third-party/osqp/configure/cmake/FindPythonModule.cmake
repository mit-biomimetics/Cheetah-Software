# Find python module and version
#
# It sets the variables (given module called MODULE)
# unless MODULE is __FUTURE__. In that case
# it only checks if it has been found.
#
# MODULE_FOUND               - has the module been found?
# MODULE_VERSION             - module version as a string
# MODULE_VERSION_MAJOR       - major version number
# MODULE_VERSION_MINOR       - minor version number
# MODULE_VERSION_PATCH       - patch version number
# MODULE_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
#


function(find_python_module module)
    # Write module in upper and lower case
	string(TOUPPER ${module} module_upper)
	string(TOLOWER ${module} module_lower)

    unset(${module_upper}_VERSION)
    #
    # if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
    #     set(${module_upper}_FIND_REQUIRED TRUE)
    # endif()

    if(PYTHONINTERP_FOUND)
        if (NOT ${module} STREQUAL "__future__")
            execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
              "import ${module_lower} as n; print(n.__version__);"
              RESULT_VARIABLE __result
              OUTPUT_VARIABLE __output
              ERROR_QUIET
              OUTPUT_STRIP_TRAILING_WHITESPACE)

              if(__result MATCHES 0)
                string(REGEX REPLACE ";" "\\\\;" __values ${__output})
                string(REGEX REPLACE "\r?\n" ";"    __values ${__values})
                list(GET __values 0 ${module_upper}_VERSION)

                string(REGEX MATCH "^([0-9])+\\.([0-9])+\\.([0-9])+" __ver_check "${${module_upper}_VERSION}")

                if(NOT "${__ver_check}" STREQUAL "")
                  set(${module_upper}_VERSION_MAJOR ${CMAKE_MATCH_1})
                  set(${module_upper}_VERSION_MINOR ${CMAKE_MATCH_2})
                  set(${module_upper}_VERSION_PATCH ${CMAKE_MATCH_3})
                  math(EXPR ${module_upper}_VERSION_DECIMAL
                    "(${CMAKE_MATCH_1} * 10000) + (${CMAKE_MATCH_2} * 100) + ${CMAKE_MATCH_3}")
                else()
                 unset(${module_upper}_VERSION)
                 message(STATUS "Requested ${module_lower} version, but got instead:\n${__output}\n")
                endif()

                find_package_handle_standard_args(${module_upper}
                    FOUND_VAR ${module_upper}_FOUND
                    REQUIRED_VARS ${module_upper}_VERSION
                    VERSION_VAR  ${module_upper}_VERSION)

              endif()
        else()
            execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
              "import ${module_lower} as n"
              RESULT_VARIABLE __result
              OUTPUT_VARIABLE __output
              ERROR_QUIET
              OUTPUT_STRIP_TRAILING_WHITESPACE)
              if(NOT __result)
                  set(${module_upper}_FOUND ON)
              endif()
              message(STATUS "Found Python __FUTURE__")
        endif()


    else()
        message(STATUS "Python interpreter not found. To find ${module} you need the Python interpreter.")
    endif()


    # Set variables in parent scope
    if(${module_upper}_FOUND)
        set(${module_upper}_FOUND ON PARENT_SCOPE)
        if (NOT ${module} STREQUAL "__future__")
            set(${module_upper}_VERSION ${${module_upper}_VERSION} PARENT_SCOPE)
            set(${module_upper}_VERSION_MAJOR ${${module_upper}_VERSION_MAJOR} PARENT_SCOPE)
            set(${module_upper}_VERSION_MINOR ${${module_upper}_VERSION_MINOR} PARENT_SCOPE)
            set(${module_upper}_VERSION_PATCH ${${module_upper}_VERSION_PATCH} PARENT_SCOPE)
            set(${module_upper}_VERSION_DECIMAL ${${module_upper}_VERSION_DECIMAL} PARENT_SCOPE)
        endif()
    endif()

    # Clear variables
    osqp_clear_vars(__result __output __values __ver_check)

endfunction(find_python_module)
