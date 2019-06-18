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


#ifndef CASADI_SNOPT_INTERFACE_HPP
#define CASADI_SNOPT_INTERFACE_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/interfaces/snopt/casadi_nlpsol_snopt_export.h"
extern "C" {
#include "snopt_cwrap.h" // NOLINT(build/include)
}

/** \defgroup plugin_Nlpsol_snopt
  SNOPT interface
*/

/** \pluginsection{Nlpsol,snopt} */

/// \cond INTERNAL
namespace casadi {

  // Forward declaration
  class SnoptInterface;

  struct CASADI_NLPSOL_SNOPT_EXPORT SnoptMemory : public NlpsolMemory {
    /// Function object
    const SnoptInterface& self;

    // Current solution
    double *xk2, *lam_gk, *lam_xk;

    // Current calculated quantities
    double *gk, *jac_fk, *jac_gk;

    std::vector<double> bl, bu, xx;

    std::vector<int> hs, locJ, indJ;

    casadi_int n_iter; // number of major iterations

    std::vector<double> A_data, valJ, rc, pi;

    // Memory pool
    static std::vector<SnoptMemory*> mempool;
    int memind;

    int return_status;

    /// Constructor
    SnoptMemory(const SnoptInterface& self);

    /// Destructor
    ~SnoptMemory();
  };

  /** \brief \pluginbrief{Nlpsol,snopt}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_snopt
  */
  class CASADI_NLPSOL_SNOPT_EXPORT SnoptInterface : public Nlpsol {
  public:
    // NLP functions
    Function f_fcn_;
    Function g_fcn_;
    Function jac_g_fcn_;
    Function jac_f_fcn_;
    Function gf_jg_fcn_;
    Function hess_l_fcn_;
    Sparsity jacg_sp_;

    // Constructor
    explicit SnoptInterface(const std::string& name, const Function& nlp);

    // Destructor
    ~SnoptInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "snopt";}

    // Get name of the class
    std::string class_name() const override { return "SnoptInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new SnoptInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SnoptMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SnoptMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /// Exact Hessian?
    bool exact_hessian_;

    static std::map<int, std::string> status_;
    static std::map<int, std::string> secondary_status_;

    std::string formatStatus(int status) const;
    std::string formatSecondaryStatus(int status) const;

    void userfun(SnoptMemory* m, int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
                 double* x, double* fObj, double*gObj, double* fCon, double* gCon,
                 int nState, char* cu, int lencu, int* iu, int leniu, double* ru, int lenru) const;

    casadi_int nnJac_;
    casadi_int nnObj_;
    casadi_int nnCon_;

    IM A_structure_;

    casadi_int m_;
    casadi_int iObj_;

    static void userfunPtr(int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac,
                           double *x, double *fObj, double *gObj, double * fCon, double* gCon,
                           int* nState, char* cu, int* lencu, int* iu, int* leniu,
                           double* ru, int *lenru);

    // Matrix A has a linear objective row
    bool jacF_row_;
    // Matrix A has a dummy row
    bool dummyrow_;

    /// A documentation string
    static const std::string meta_doc;

    /// Warm-start settings
    casadi_int Cold_;

    double inf_;

  private:
      // options
      Dict opts_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SNOPT_INTERFACE_HPP
