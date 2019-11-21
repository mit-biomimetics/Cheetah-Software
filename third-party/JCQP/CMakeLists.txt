include_directories(amd/include)
add_library(JCQP SHARED
        QpProblem.cpp
        ProblemGenerator.cpp
        eigenvalues.cpp
        CholeskyDenseSolver.cpp
        CholeskySparseSolver.cpp
        SparseMatrixMath.cpp
        #OsqpTest.cpp
        #qpOASES.cpp
        #ThreadGroup.cpp
        amd/src/amd_1.c
        amd/src/amd_2.c
        amd/src/amd_aat.c
        amd/src/amd_control.c
        amd/src/amd_defaults.c
        amd/src/amd_info.c
        amd/src/amd_order.c
        amd/src/amd_post_tree.c
        amd/src/amd_postorder.c
        amd/src/amd_preprocess.c
        amd/src/amd_valid.c
        amd/src/SuiteSparse_config.c)


#target_link_libraries(JCQP OsqpEigen::OsqpEigen osqp::osqp qpOASES pthread)
target_link_libraries(JCQP pthread)
