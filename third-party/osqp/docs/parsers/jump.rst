JuMP
=====

`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ supports the OSQP solver using the `MathProgBase interface <https://github.com/JuliaOpt/MathProgBase.jl>`_. 
You can define a JuMP model to be solved via OSQP as follows


.. code:: julia

   # Load JuMP and OSQP
   using JuMP, OSQP

   # Create OSQP Solver instance
   s = OSQPMathProgBaseInterface.OSQPSolver(verbose=false)

   # Create JuMP model
   model = Model(solver=s)

   ...

Note that here we set the verbosity to :code:`false`.
After defining your model, you can solve it by just calling

.. code:: julia

   solve(model)


For more details on how to create and modify the models, see the `JuMP Documentation <https://github.com/JuliaOpt/JuMP.jl>`_.