if(WITH_TINYXML AND WITH_BUILD_TINYXML)
  add_subdirectory(tinyxml-2.6.2)
endif()

if(WITH_CSPARSE AND WITH_BUILD_CSPARSE)
  add_subdirectory(CSparse)
endif()

if(WITH_SUNDIALS AND WITH_BUILD_SUNDIALS)
  add_subdirectory(casadi-sundials)
endif()

if(WITH_SNOPT)
  add_subdirectory(snopt-interface)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpOASES)
endif()

if(WITH_BUILD_DSDP AND DSDP_FOUND)
  add_subdirectory(DSDP)
endif()

if(WITH_HPMPC AND WITH_BUILD_HPMPC)
  add_subdirectory(hpmpc)
endif()

if(WITH_BLASFEO AND WITH_BUILD_BLASFEO)
  add_subdirectory(blasfeo)
endif()
