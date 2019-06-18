Contributing
=============

OSQP is an open-source project open to any academic or commercial applications.
Contributions are welcome as `GitHub pull requests <https://help.github.com/articles/creating-a-pull-request/>`_ in any part of the project such as

* algorithm developments
* interfaces to other languages
* compatibility with new architectures
* linear system solvers interfaces


.. _interfacing_new_linear_system_solvers :

Interfacing new linear system solvers
-------------------------------------
OSQP is designed to be easily interfaced to new linear system solvers via dynamic library loading.
To add a linear system solver interface you need to edit the :code:`lin_sys/` directory subfolder :code:`direct/` or :code:`indirect/` depending on the type of solver.
Create a subdirectory with your solver name with four files:

* Dynamic library loading: :code:`mysolver_loader.c` and :code:`mysolver_loader.h`.
* Linear system solution: :code:`mysolver.c` and :code:`mysolver.h`.

We suggest you to have a look at the `MKL Pardiso solver interface <https://github.com/oxfordcontrol/osqp/tree/master/lin_sys/direct/pardiso>`_ for more details.

Dynamic library loading
^^^^^^^^^^^^^^^^^^^^^^^
In this part define the methods to load the shared library and obtain the functions required to solve the linear system.
The main functions to be exported are :code:`lh_load_mysolver(const char* libname)` and :code:`lh_unload_mysolver()`.
In addition, the file :code:`mysolver_loader.c` must define static function pointers to the shared library functions to be loaded.

Linear system solution
^^^^^^^^^^^^^^^^^^^^^^
In this part we define the core of the interface: **linear system solver object**.
The main functions are the external method

* :code:`init_linsys_solver_mysolver`: create the instance and perform the setup

and the internal methods of the object

* :code:`free_linsys_solver_mysolver`: free the instance
* :code:`solve_linsys_mysolver`: solve the linear system
* :code:`update_matrices`: update problem matrices
* :code:`update_rho_vec`: update :math:`\rho` as a diagonal vector.

After the initializations these functions are assigned to the internal pointers so that, for an instance :code:`s` they can be called as :code:`s->free`, :code:`s->solve`, :code:`s->update_matrices` and :code:`s->update_rho_vec`.

The linear system solver object is defined in :code:`mysolver.h` as follows

.. code:: c

        typedef struct mysolver mysolver_solver;

        struct mysolver {
            // Methods
            enum linsys_solver_type type; // Linear system solver defined in constants.h

            c_int (*solve)(struct mysolver * self, c_float * b, const OSQPSettings * settings);
            void (*free)(struct mysolver * self);
            c_int (*update_matrices)(struct mysolver * self, const csc *P, const csc *A, const OSQPSettings *settings);
            c_int (*update_rho_vec)(struct mysolver * self, const c_float * rho_vec, const c_int m);

            // Attributes
            c_int nthreads; // Number of threads used (required!)

            // Internal attributes of the solver
            ...

            // Internal attributes required for matrix updates
            c_int * Pdiag_idx, Pdiag_n;  ///< index and number of diagonal elements in P
            c_int * PtoKKT, * AtoKKT;    ///< Index of elements from P and A to KKT matrix
            c_int * rhotoKKT;            ///< Index of rho places in KKT matrix
            ...

        };

        // Initialize mysolver solver
        mysolver_solver *init_linsys_solver_mysolver(const csc * P, const csc * A, c_float sigma, c_float * rho_vec, c_int polish);

        // Solve linear system and store result in b
        c_int solve_linsys_mysolver(mysolver_solver * s, c_float * b, const OSQPSettings * settings);

         // Update linear system solver matrices
        c_int update_linsys_solver_matrices_mysolver(mysolver_solver * s,
                        const csc *P, const csc *A, const OSQPSettings *settings);

        // Update rho parameter in linear system solver structure
        c_int update_linsys_solver_rho_vec_mysolver(mysolver_solver * s, const c_float * rho_vec, const c_int m);

        // Free linear system solver
        void free_linsys_solver_mysolver(mysolver_solver * s);


The function details are coded in the :code:`mysolver.c` file.
