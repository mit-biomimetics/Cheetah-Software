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


#ifndef CASADI_IMPLICIT_TO_NLP_HPP
#define CASADI_IMPLICIT_TO_NLP_HPP

#include "casadi/core/rootfinder_impl.hpp"
#include <casadi/solvers/casadi_rootfinder_nlpsol_export.h>

/** \defgroup plugin_Rootfinder_nlp
  Use an Nlpsol as Rootfinder solver
*/
/** \pluginsection{Rootfinder,nlpsol} */

/// \cond INTERNAL
namespace casadi {
  // Memory
  struct CASADI_ROOTFINDER_NLPSOL_EXPORT ImplicitToNlpMemory
    : public RootfinderMemory {

    // Bounds
    double *lbx, *ubx;
    // Parameters
    double *p;
    // solution
    double *x;
  };

  /** \brief  \pluginbrief{Rootfinder,nlp}

   @copydoc Rootfinder_doc
   @copydoc plugin_Rootfinder_nlp

   \author Joris Gillis
   \date 2012
  */
  class CASADI_ROOTFINDER_NLPSOL_EXPORT ImplicitToNlp : public Rootfinder {
  public:
    /** \brief  Constructor */
    explicit ImplicitToNlp(const std::string& name, const Function& f);

    /** \brief  Destructor */
    ~ImplicitToNlp() override;

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new ImplicitToNlp(name, f);
    }

    // Get name of the plugin
    const char* plugin_name() const override { return "nlpsol";}

    // Name of the class
    std::string class_name() const override { return "ImplicitToNlp";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new ImplicitToNlpMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<ImplicitToNlpMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// Solve the system of equations and calculate derivatives
    int solve(void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /// NLP solver
    Function solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_IMPLICIT_TO_NLP_HPP
