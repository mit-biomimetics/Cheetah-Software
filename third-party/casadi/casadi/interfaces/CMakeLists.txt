if(WITH_LAPACK)
  add_subdirectory(lapack)
endif()

if(WITH_SUNDIALS)
  add_subdirectory(sundials)
endif()

if(WITH_IPOPT)
  add_subdirectory(ipopt)
endif()

if(WITH_BONMIN)
  add_subdirectory(bonmin)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpoases)
endif()

if(WITH_KNITRO)
  add_subdirectory(knitro)
endif()

if(WITH_CPLEX)
  add_subdirectory(cplex)
endif()

if(WITH_CLP)
  add_subdirectory(clp)
endif()

if(WITH_CSPARSE)
  add_subdirectory(csparse)
endif()

if(WITH_CLANG)
  add_subdirectory(clang)
endif()

if(WITH_HSL)
  add_subdirectory(hsl)
endif()

if(WITH_OOQP)
  add_subdirectory(ooqp)
endif()

if(WITH_GUROBI)
  add_subdirectory(gurobi)
endif()

if(WITH_WORHP)
  add_subdirectory(worhp)
endif()

if(WITH_SQIC)
  add_subdirectory(sqic)
endif()

if(WITH_SNOPT)
  add_subdirectory(snopt)
endif()

if(WITH_TINYXML)
  add_subdirectory(tinyxml)
endif()

if(WITH_BLOCKSQP)
  add_subdirectory(blocksqp)
endif()

if(WITH_HPMPC)
  add_subdirectory(hpmpc)
endif()

if(WITH_SLICOT)
  add_subdirectory(slicot)
endif()

if(WITH_AMPL)
  add_subdirectory(ampl)
endif()

if(WITH_LINT)
  set(LINT_TARGETS ${LINT_TARGETS} PARENT_SCOPE)
endif()

if(WITH_SPELL)
  set(SPELL_TARGETS ${SPELL_TARGETS} PARENT_SCOPE)
endif()

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
