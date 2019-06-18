set(test_headers
	${CMAKE_CURRENT_SOURCE_DIR}/test_basic.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_identity.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_rank_deficient.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_singleton.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_sym_structure.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_tril_structure.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_two_by_two.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_zero_on_diag.h
	${CMAKE_CURRENT_SOURCE_DIR}/test_osqp_kkt.h
	PARENT_SCOPE)

# Include this directory for test headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
