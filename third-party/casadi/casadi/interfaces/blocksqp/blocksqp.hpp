/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_BLOCKSQP_HPP
#define CASADI_BLOCKSQP_HPP

#include <casadi/interfaces/blocksqp/casadi_nlpsol_blocksqp_export.h>
#include "casadi/core/linsol.hpp"
#include "casadi/core/nlpsol_impl.hpp"
#include "../qpoases/qpoases_interface.hpp"

/** \defgroup plugin_Nlpsol_blocksqp
  * This is a modified version of blockSQP by Janka et al.
  *
  * \author Dennis Janka, Joel Andersson
  * \date 2012-2015, 2016
*/

/** \pluginsection{Nlpsol,blocksqp} */

/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class Blocksqp;

  struct CASADI_NLPSOL_BLOCKSQP_EXPORT BlocksqpMemory : public NlpsolMemory {
    /// Constructor
    BlocksqpMemory();

    /// Destructor
    ~BlocksqpMemory();

    qpOASES::SymSparseMat *H;
    qpOASES::Matrix *A;
    qpOASES::SQProblem* qp;

    // [Workaround] qpOASES memory block
    QpoasesMemory* qpoases_mem;

    // Stats
    casadi_int itCount;  // iteration number
    casadi_int qpIterations;  // number of qp iterations in the current major iteration
    casadi_int qpIterations2;  // number of qp iterations for solving convexified QPs
    casadi_int qpItTotal;  // total number of qp iterations
    casadi_int qpResolve;  // how often has QP to be convexified and resolved?
    casadi_int nFunCalls;  // number of function calls
    casadi_int nDerCalls;  // number of derivative calls
    casadi_int nRestHeurCalls;  // number calls to feasibility restoration heuristic
    casadi_int nRestPhaseCalls;  // number calls to feasibility restoration phase
    casadi_int rejectedSR1;  // count how often the SR1 update is rejected
    casadi_int hessSkipped;  // number of block updates skipped in the current iteration
    casadi_int hessDamped;  // number of block updates damped in the current iteration
    casadi_int nTotalUpdates;
    casadi_int nTotalSkippedUpdates;
    double averageSizingFactor;  // average value (over all blocks) of COL sizing factor

    // Variables that are updated during one SQP iteration
    double obj;  // objective value
    double qpObj;  // objective value of last QP subproblem
    double cNorm;  // constraint violation
    double cNormS;  // scaled constraint violation
    double gradNorm;  // norm of Lagrangian gradient
    double lambdaStepNorm;  // norm of step in dual variables
    double tol;  // current optimality tolerance

    double *lam_xk, *lam_gk; // dual variables
    double *gk;  // constraint vector

    double *jac_g;  // nonzero elements of Jacobian (length)

    double *deltaMat;  // last m primal steps
    double *dxk;  // alias for current step
    double *grad_fk;  // gradient of objective
    double *grad_lagk;  // gradient of Lagrangian
    double *gammaMat;  // Lagrangian gradient differences for last m steps
    double *gamma;  // alias for current Lagrangian gradient

    double **hess;  // [blockwise] pointer to current Hessian of the Lagrangian
    double **hess1;  // [blockwise] first Hessian approximation
    double **hess2;  // [blockwise] second Hessian approximation (convexified)
    double *hess_lag;  // nonzero elements of Hessian (length)
    int *hessIndRow;  // row indices (length)
    int *hessIndCol;  // indices to first entry of columns (nCols+1)
    int *hessIndLo;  // Indices to first entry of lower triangle (including diagonal) (nCols)

    /*
     * Variables for QP solver
     */
    double *lbx_qp, *ubx_qp, *lba_qp, *uba_qp;  // bounds for current step
    double* lam_qp;  // dual variables of QP
    double* jac_times_dxk;  // product of constraint Jacobian with dxk

    /*
     * For modified BFGS updates
     */
    double *delta_norm;  // sTs
    double* delta_norm_old;  // (from previous iteration)
    double* delta_gamma;  // sTy
    double* delta_gamma_old;  // (from previous iteration)
    casadi_int *noUpdateCounter;  // count skipped updates for each block

    /*
     * Variables for globalization strategy
     */
    casadi_int steptype;  // is current step a restoration step (1)?
    double alpha;  // stepsize for line search
    casadi_int nSOCS;  // number of second-order correction steps
    casadi_int reducedStepCount;  // count number of consecutive reduced steps,
    double* delta_h; // inertia correction (filter line search w indef Hessian)
    double* trial_xk;  // new trial iterate (for line search)
    std::set< std::pair<double, double> > filter; // Filter contains pairs (constrVio, objective)

    std::vector<int> colind, row;

    // Temporary memory
    double* jac;
  };

  /** \brief \pluginbrief{Nlpsol,blocksqp}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_blocksqp
  */
  class CASADI_NLPSOL_BLOCKSQP_EXPORT Blocksqp : public Nlpsol {
  public:
    explicit Blocksqp(const std::string& name, const Function& nlp);
    ~Blocksqp() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "blocksqp";}

    // Get name of the class
    std::string class_name() const override { return "Blocksqp";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Blocksqp(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new BlocksqpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<BlocksqpMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Block partitioning
    casadi_int nblocks_;
    std::vector<casadi_int> blocks_;
    std::vector<casadi_int> dim_;
    casadi_int nnz_H_;

    // Jacobian/Hessian sparsity
    Sparsity Asp_, Hsp_;

    /// Main Loop of SQP method
    casadi_int run(BlocksqpMemory* m, casadi_int maxIt, casadi_int warmStart = 0) const;
    /// Compute gradient of Lagrangian function (sparse version)
    void calcLagrangeGradient(BlocksqpMemory* m,
      const double* lam_x, const double* lam_g,
      const double* grad_f, double *jacNz,
      double *grad_lag, casadi_int flag) const;

    /// Overloaded function for convenience, uses current variables of SQPiterate vars
    void calcLagrangeGradient(BlocksqpMemory* m, double* grad_lag, casadi_int flag) const;
    /// Print information about the SQP method
    void printInfo(BlocksqpMemory* m) const;
    /// Update optimization tolerance (similar to SNOPT) in current iterate
    bool calcOptTol(BlocksqpMemory* m) const;

    /*
     * Solve QP subproblem
     */
    // Update the bounds on the current step, i.e. the QP variables
    void updateStepBounds(BlocksqpMemory* m, bool soc) const;
    // Solve a QP with QPOPT or qpOASES to obtain a step deltaXi and estimates
    // for the Lagrange multipliers
    casadi_int solveQP(BlocksqpMemory* m, double* deltaXi, double* lambdaQP,
      bool matricesChanged = true) const;
    // Compute the next Hessian in the inner loop of increasingly convexified
    // QPs and store it in vars->hess2
    void computeNextHessian(BlocksqpMemory* m, casadi_int idx, casadi_int maxQP) const;

    /*
     * Globalization Strategy
     */
    /// No globalization strategy
    casadi_int fullstep(BlocksqpMemory* m) const;
    /// Set new primal dual iterate
    void acceptStep(BlocksqpMemory* m, const double* deltaXi,
      const double* lambdaQP, double alpha, casadi_int nSOCS) const;
    // Overloaded function for convenience, uses current variables of SQPiterate vars
    void acceptStep(BlocksqpMemory* m, double alpha) const;
    // Reduce stepsize if a step is rejected
    void reduceStepsize(BlocksqpMemory* m, double *alpha) const;
    // Determine steplength alpha by a filter based line search similar to IPOPT
    casadi_int filterLineSearch(BlocksqpMemory* m) const;
    // Remove all entries from filter
    void initializeFilter(BlocksqpMemory* m) const;
    // Is a pair (cNorm, obj) in the current filter?
    bool pairInFilter(BlocksqpMemory* m, double cNorm, double obj) const;
    // Augment current filter by pair (cNorm, obj)
    void augmentFilter(BlocksqpMemory* m, double cNorm, double obj) const;
    // Perform a second order correction step (solve QP)
    bool secondOrderCorrection(BlocksqpMemory* m, double cNorm, double cNormTrial,
      double dfTdeltaXi, bool swCond, casadi_int it) const;
    // Reduce stepsize if a second order correction step is rejected
    void reduceSOCStepsize(BlocksqpMemory* m, double *alphaSOC) const;
    // Start feasibility restoration heuristic
    casadi_int feasibilityRestorationHeuristic(BlocksqpMemory* m) const;
    // Start feasibility restoration phase (solve NLP)
    casadi_int feasibilityRestorationPhase(BlocksqpMemory* m) const;
    // Check if full step reduces KKT error
    casadi_int kktErrorReduction(BlocksqpMemory* m) const;

    /*
     * Hessian Approximation
     */
    // Set initial Hessian: Identity matrix
    void calcInitialHessian(BlocksqpMemory* m) const;
    // [blockwise] Set initial Hessian: Identity matrix
    void calcInitialHessian(BlocksqpMemory* m, casadi_int b) const;
    // Reset Hessian to identity and remove past information on Lagrange gradient and steps
    void resetHessian(BlocksqpMemory* m) const;
    // [blockwise] Reset Hessian to identity and remove past information on
    // Lagrange gradient and steps
    void resetHessian(BlocksqpMemory* m, casadi_int b) const;
    // Compute full memory Hessian approximations based on update formulas
    void calcHessianUpdate(BlocksqpMemory* m, casadi_int updateType, casadi_int hessScaling) const;
    // Compute limited memory Hessian approximations based on update formulas
    void calcHessianUpdateLimitedMemory(BlocksqpMemory* m,
        casadi_int updateType, casadi_int hessScaling) const;
    // [blockwise] Compute new approximation for Hessian by SR1 update
    void calcSR1(BlocksqpMemory* m, const double* gamma, const double* delta,
      casadi_int b) const;
    // [blockwise] Compute new approximation for Hessian by BFGS update with Powell modification
    void calcBFGS(BlocksqpMemory* m, const double* gamma, const double* delta,
      casadi_int b) const;
    // Set pointer to correct step and Lagrange gradient difference in a limited memory context
    void updateDeltaGamma(BlocksqpMemory* m) const;

    /*
     * Scaling of Hessian Approximation
     */
    // [blockwise] Size Hessian using SP, OL, or mean sizing factor
    void sizeInitialHessian(BlocksqpMemory* m, const double* gamma,
      const double* delta, casadi_int b, casadi_int option) const;
    // [blockwise] Size Hessian using the COL scaling factor
    void sizeHessianCOL(BlocksqpMemory* m, const double* gamma,
      const double* delta, casadi_int b) const;

    /*
    * STATS
    */
    void initStats(BlocksqpMemory* m) const;
    void updateStats(BlocksqpMemory* m) const;
    /// Print one line of output to stdout about the current iteration
    void printProgress(BlocksqpMemory* m) const;
    /// Reset variables that any SQP code needs
    void reset_sqp(BlocksqpMemory* m) const;
    /// Convert *hess to column compressed sparse format
    void convertHessian(BlocksqpMemory* m) const;
    /// Set initial filter, objective function, tolerances etc.
    void initIterate(BlocksqpMemory* m) const;

    /// Evaluate objective and constraints, including derivatives
    casadi_int evaluate(BlocksqpMemory* m, double *f, double *g,
                 double *grad_f, double *jac_g) const;

    /// Evaluate objective and constraints, no derivatives
    casadi_int evaluate(BlocksqpMemory* m, const double *xk,
                 double *f, double *g) const;

    //  Declaration of general purpose routines for matrix and vector computations
    double lInfConstraintNorm(BlocksqpMemory* m, const double* xk, const double* g) const;

    /// QP solver for the subproblems
    //Function qpsol_;

    // Linear solver in qpOASES
    std::string linsol_plugin_;

    // General options
    bool print_header_;
    bool print_iteration_;
    double eps_;  // values smaller than this are regarded as numerically zero
    double opttol_;  // optimality tolerance
    double nlinfeastol_; // nonlinear feasibility tolerance

    // Algorithmic options
    bool schur_;  // Use qpOASES schur compliment approach
    bool globalization_; // Globalization strategy
    bool restore_feas_;// Use feasibility restoration phase
    casadi_int max_line_search_;  // Maximum number of steps in line search
    casadi_int max_consec_reduced_steps_;// Maximum number of consecutive reduced steps
    casadi_int max_consec_skipped_updates_; // Maximum number of consecutive skipped updates
    casadi_int max_it_qp_;  // Maximum number of QP iterations per SQP iteration
    casadi_int max_iter_; // Maximum number of SQP steps
    bool warmstart_; // Use warmstarting
    bool block_hess_;  // Blockwise Hessian approximation?
    casadi_int hess_scaling_;// Scaling strategy for Hessian approximation
    casadi_int fallback_scaling_;  // If indefinite update is used, the type of fallback strategy
    double max_time_qp_;  // Maximum number of time in seconds per QP solve per SQP iteration
    double ini_hess_diag_;  // Initial Hessian guess: diagonal matrix diag(iniHessDiag)
    double col_eps_;  // epsilon for COL scaling strategy
    double col_tau1_; // tau1 for COL scaling strategy
    double col_tau2_; // tau2 for COL scaling strategy
    casadi_int hess_damp_;  // activate Powell damping for BFGS
    double hess_damp_fac_;  // damping factor for BFGS Powell modification
    casadi_int hess_update_; // Type of Hessian approximation
    casadi_int fallback_update_;  // If indefinite update is used, the type of fallback strategy
    casadi_int hess_lim_mem_; // Full or limited memory
    casadi_int hess_memsize_;  // Memory size for L-BFGS updates
    casadi_int which_second_derv_; // For which block should second derivatives be provided
                                   // by the user
    bool skip_first_globalization_;  // No globalization strategy in first iteration
    casadi_int conv_strategy_;  // Convexification strategy
    casadi_int max_conv_qp_;  // How many additional QPs may be solved for convexification
                              // per iteration?

    // Filter line search parameters, cf. IPOPT paper
    casadi_int max_soc_iter_; // Maximum number of SOC line search iterations
    double gamma_theta_;
    double gamma_f_;
    double kappa_soc_;
    double kappa_f_;
    double theta_max_;
    double theta_min_;
    double delta_;
    double s_theta_;
    double s_f_;
    double kappa_minus_;
    double kappa_plus_;
    double kappa_plus_max_;
    double delta_h0_;
    double eta_;
    double obj_lo_;
    double obj_up_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BLOCKSQP_HPP
