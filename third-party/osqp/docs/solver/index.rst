The solver
===========

Problem statement
-----------------

OSQP solves convex quadratic programs (QPs) of the form

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T P x + q^T x \\
    \mbox{subject to} & l \leq A x \leq u
  \end{array}

where :math:`x\in\mathbf{R}^{n}` is the optimization variable.
The objective function is defined by a positive semidefinite matrix
:math:`P \in \mathbf{S}^{n}_{+}` and vector :math:`q\in \mathbf{R}^{n}`.
The linear constraints are defined by matrix :math:`A\in\mathbf{R}^{m \times n}`
and vectors :math:`l \in \mathbf{R}^{m} \cup \{-\infty\}^{m}`,
:math:`u \in \mathbf{R}^{m} \cup \{+\infty\}^{m}`.


Algorithm
-------------------------

The solver runs the following `ADMM algorithm <http://web.stanford.edu/~boyd/papers/admm_distr_stats.html>`_ (for more details see the related papers at the :ref:`citing` section):


.. math::

   (x^{k+1}, \nu^{k+1}) & \gets \text{solve linear system}\\
   \tilde{z}^{k+1} & \gets z^{k} + \rho^{-1}(\nu^{k+1} - y^{k})\\
   z^{k+1} &\gets \Pi(\tilde{z}^{k} + \rho^{-1}y^{k})\\
   y^{k+1} &\gets y^{k} + \rho (\tilde{z}^{k+1} - z^{k+1})

where :math:`\Pi` is the projection onto the hyperbox :math:`[l,u]`. 
:math:`\rho` is the ADMM step-size. 


Linear system solution
^^^^^^^^^^^^^^^^^^^^^^^
The linear system solution is the core part of the algorithm.
It can be done using a **direct** or **indirect** method.

With a **direct linear system solver** we solve the following linear system with a quasi-definite matrix 

.. math::

   \begin{bmatrix} P + \sigma I & A^T \\ A & -\rho^{-1}I \end{bmatrix} \begin{bmatrix} x^{k+1} \\ \nu^{k+1} \end{bmatrix}= \begin{bmatrix}\sigma x^{k} - q \\ z^{k} - \rho^{-1} y^k \end{bmatrix}.

With an **indirect linear system solver** we solve the following linear system with a positive definite matrix 


.. math::

	\left(P + \sigma I + \rho A^T A \right)x^{k+1} = \sigma x^{k} - q + A^T (\rho z^{k} - y^{k}).


OSQP core is designed to support different linear system solvers.
For their installation see :ref:`this section <linear_system_solvers_installation>`.
To specify your favorite linear system solver see :ref:`this section <linear_system_solvers_setting>`.

Convergence
^^^^^^^^^^^
At each iteration :math:`k` OSQP produces a tuple :math:`(x^{k}, z^{k}, y^{k})` with :math:`x^{k} \in \mathbf{R}^{n}` and :math:`z^{k}, y^{k} \in \mathbf{R}^{m}`.

The primal and and dual residuals associated to :math:`(x^{k}, z^{k}, y^{k})` are

.. math::

   \begin{aligned}
   r_{\rm prim}^{k} &= Ax^{k} - z^{k}\\
   r_{\rm dual}^{k} &= Px^{k} + q + A^{T} y^{k}.
   \end{aligned}

Complementary slackness is satisfied by construction at machine precision. If the problem is feasible, the residuals converge to zero as :math:`k\to\infty`. The algorithm stops when the norms of :math:`r_{\rm prim}^{k}` and :math:`r_{\rm dual}^{k}` are within the specified tolerance levels :math:`\epsilon_{\rm prim}>0` and :math:`\epsilon_{\rm dual}>0` 

.. math::

    \| r_{\rm prim}^{k} \|_{\infty} \le 
   \epsilon_{\rm prim},
    \quad
    \| r_{\rm dual}^{k} \|_{\infty} \le 
   \epsilon_{\rm dual}.

We set the tolerance levels as
    
.. math:: 

    \epsilon_{\rm prim} &= \epsilon_{\rm abs} + \epsilon_{\rm rel} \max\lbrace \|Ax^{k}\|_{\infty}, \| z^{k} \|_{\infty} \rbrace \\
    \epsilon_{\rm dual} &= \epsilon_{\rm abs} + \epsilon_{\rm rel} \max\lbrace \| P x^{k} \|_{\infty}, \| A^T y^{k} \|_{\infty}, \| q \|_{\infty} \rbrace.


.. _rho_step_size :

:math:`\rho` step-size 
^^^^^^^^^^^^^^^^^^^^^^^^^
To ensure quick convergence of the algorithm we adapt :math:`\rho` by balancing the residuals. 
In default mode, the inteval (*i.e.*, number of iterations) at which we update :math:`\rho` is defined by a time measurement.
When the iterations time becomes greater than a certain fraction of the setup time, *i.e.* :code:`adaptive_rho_fraction`, we set the current number of iterations as the interval to update :math:`\rho`.
The update happens as follows 

.. math::

    \rho^{k+1} \gets \rho^{k} \sqrt{\frac{\|r_{\rm prim}\|_{\infty}}{\|r_{\rm dual}\|_{\infty}}}.


Note that :math:`\rho` is updated only if it is sufficiently different than the current one.
In particular if it is :code:`adaptive_rho_tolerance` times larger or smaller than the current one.



Infeasible problems
-------------------------------

OSQP is able to detect if the problem is primal or dual infeasible.


Primal infeasibility
^^^^^^^^^^^^^^^^^^^^

When the problem is primal infeasible, the algorithm generates a certificate of infeasibility, *i.e.*, a vector :math:`v\in\mathbf{R}^{m}` such that

.. math::

   A^T v = 0, \quad u^T v_{+} + l^T v_{-} < 0,

where :math:`v_+=\max(v,0)` and :math:`v_-=\min(v,0)`.

The primal infeasibility check is then

.. math::

	\left\|A^T v \right\|_{\infty} \le \epsilon_{\rm prim\_inf}, \quad u^T (v)_{+} + l^T (v)_{-} \le - \epsilon_{\rm prim\_inf}.



Dual infeasibility
^^^^^^^^^^^^^^^^^^^^

When the problem is dual infeasible, OSQP generates a vector :math:`s\in\mathbf{R}^{n}` being a certificate of dual infeasibility, *i.e.*,

.. math::

   P s = 0, \quad q^T s < 0, \quad (As)_i = \begin{cases} 0 & l_i \in \mathbf{R}, u_i\in\mathbf{R} \\ \ge 0 & l_i\in\mathbf{R}, u_i=+\infty \\ \le 0 & u_i\in\mathbf{R}, l_i=-\infty. \end{cases}


The dual infeasibility check is then

.. math::

        \| P s \|_{\infty} \le \epsilon_{\rm dual\_inf} , \quad
        q^T s \le -\epsilon_{\rm dual\_inf}, \\
        (A s)_i \begin{cases} \in \left[-\epsilon_{\rm dual\_inf}, \epsilon_{\rm dual\_inf}\right] & u_i, l_i \in \mathbf{R}\\
        \ge \epsilon_{\rm dual\_inf} &u_i = +\infty\\
        \le  -\epsilon_{\rm dual\_inf} &l_i = -\infty.\end{cases}




Polishing
^^^^^^^^^^^

Polishing is an additional algorithm step where OSQP tries to compute a high-accuracy solution. 
We can enable it by turning the setting :code:`polish` to 1.

Polishing works by guessing the active constraints at the optimum and solving an additional linear system.
If the guess is correct, OSQP returns a high accuracy solution.
Otherwise OSQP returns the ADMM solution.
The status of the polishing phase appears in the information :code:`status_polish`.

Note that polishing requires the solution of an additional linear system and thereby, an additional factorization if the linear system solver is direct.
However, the linear system is usually much smaller than the one solved during the ADMM iterations.

The chances to have a successful polishing increase if the tolerances :code:`eps_abs` and :code:`eps_rel` are small. 
However, low tolerances might require a very large number of iterations.


