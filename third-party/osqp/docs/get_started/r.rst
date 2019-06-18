R
======

Binaries
--------

A pre-compiled version of the OSQP R interface can be installed directly from within R.   Note that this will install the OSQP interface from the current CRAN repository, which may not be the most up-to-date version:

.. code:: r

  install.packages("osqp")



The pre-compiled binaries can also be downloaded `directly from the CRAN server
<https://cran.r-project.org/web/packages/osqp/>`_.

From Sources
------------

If you would like to use the most recent version of OSQP-R and have access to git on your machine along with a suitable compiler, then you can do the following from within a terminal:

.. code:: r

  git clone --recursive https://github.com/oxfordcontrol/osqp-r.git
  cd osqp-r
  R CMD install .

From Sources (within R)
-----------------------

If you would like to install the latest version directly from with R (e.g.\ because you do not have ``git`` installed) and have a suitable compiler, then you can do the following from within R:

.. code:: r

  install.packages("remotes")
  remotes::install_github("r-lib/remotes#103")
  remotes::install_git("git://github.com/OxfordControl/osqp-r",submodules = TRUE)

Note that the second line above is necessary because the "remotes" package in R does not currently support recursive cloning of git submodules.
