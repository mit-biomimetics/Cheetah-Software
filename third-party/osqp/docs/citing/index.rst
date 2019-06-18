.. _citing :

Citing OSQP
===========

If you use OSQP for published work, we encourage you to put a star on `GitHub <https://github.com/oxfordcontrol/osqp>`_ and cite the accompanying papers:


.. glossary::

    Main paper
        Main algorithm description, derivation and benchmark available in this `preprint <https://arxiv.org/pdf/1711.08013.pdf>`_.

        .. code:: latex

          @article{osqp,
            author = {Stellato, B. and Banjac, G. and Goulart, P. and Bemporad, A. and Boyd, S.},
            title = {{OSQP}: An Operator Splitting Solver for Quadratic Programs},
            journal = {ArXiv e-prints},
            year = {2017},
            month = nov,
            adsnote = {Provided by the SAO/NASA Astrophysics Data System},
            adsurl = {http://adsabs.harvard.edu/abs/2017arXiv171108013S},
            archiveprefix = {arXiv},
            eprint = {1711.08013},
            keywords = {Mathematics - Optimization and Control},
            primaryclass = {math.OC},
          }

    Infeasibility detection
        Infeasibility detection proofs using ADMM (also for general conic programs) in this `preprint <https://stanford.edu/~boyd/papers/pdf/admm_infeas.pdf>`_.

        .. code:: latex

          @article{osqp-infeasibility,
            title   = {Infeasibility detection in the alternating direction method of multipliers for convex optimization},
            author  = {Banjac, G. and Goulart, P. and Stellato, B. and Boyd, S.},
            year    = {2017},
            url     = {https://stanford.edu/~boyd/papers/pdf/admm_infeas.pdf},
          }

    Code generation
        Code generation functionality and example in this `paper <http://stanford.edu/~boyd/papers/pdf/osqp_embedded.pdf>`_.

        .. code:: latex

          @inproceedings{osqp-codegen,
            author = {Banjac, G. and Stellato, B. and Moehle, N. and Goulart, P. and Bemporad, A. and Boyd, S.},
            title = {Embedded code generation using the {OSQP} solver},
            booktitle = {{IEEE} Conference on Decision and Control ({CDC})},
            year = {2017},
            month = dec,
            doi = {10.1109/CDC.2017.8263928},
            url = {https://doi.org/10.1109/CDC.2017.8263928},
          }
