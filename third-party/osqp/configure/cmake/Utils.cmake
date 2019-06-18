# Clears variables from list
# Usage:
#   osqp_clear_vars(<variables_list>)
macro(osqp_clear_vars)
  foreach(_var ${ARGN})
    unset(${_var})
  endforeach()
endmacro()
