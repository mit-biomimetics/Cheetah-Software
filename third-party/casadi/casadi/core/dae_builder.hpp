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


#ifndef CASADI_DAE_BUILDER_HPP
#define CASADI_DAE_BUILDER_HPP

#include "variable.hpp"

namespace casadi {

  // Forward declarations
  class XmlNode;

  /** \brief An initial-value problem in differential-algebraic equations
      <H3>Independent variables:  </H3>
      \verbatim
      t:      time
      \endverbatim

      <H3>Time-continuous variables:  </H3>
      \verbatim
      x:      states defined by ODE
      s:      implicitly defined states
      z:      algebraic variables
      u:      control signals
      q:      quadrature states
      y:      outputs
      \endverbatim

      <H3>Time-constant variables:  </H3>
      \verbatim
      p:      free parameters
      d:      dependent parameters
      \endverbatim

      <H3>Dynamic constraints (imposed everywhere):  </H3>
      \verbatim
      ODE                    \dot{x} ==  ode(t, x, s, z, u, p, d)
      DAE or implicit ODE:         0 ==  dae(t, x, s, z, u, p, d, sdot)
      algebraic equations:         0 ==  alg(t, x, s, z, u, p, d)
      quadrature equations:  \dot{q} == quad(t, x, s, z, u, p, d)
      dependent parameters:        d == ddef(t, x, s, z, u, p, d)
      output equations:            y == ydef(t, x, s, z, u, p, d)
      \endverbatim

      <H3>Point constraints (imposed pointwise):  </H3>
      \verbatim
      Initial equations:           0 == init(t, x, s, z, u, p, d, sdot)
      \endverbatim

      \date 2012-2015
      \author Joel Andersson
  */
  class CASADI_EXPORT DaeBuilder
    : public SWIG_IF_ELSE(PrintableCommon, Printable<DaeBuilder>) {
  public:

    /// Default constructor
    DaeBuilder();

    /** @name Variables and equations
     *  Public data members
     */
    ///@{
    /** \brief Independent variable (usually time) */
    MX t;

    /** \brief Differential states defined by ordinary differential equations (ODE)
     */
    std::vector<MX> x, ode, lam_ode;

    /** \brief Differential-algebraic equation (DAE) with corresponding state vector,
     * state derivatives.
     */
    std::vector<MX> s, sdot, dae, lam_dae;

    /** \brief Algebraic equations and corresponding algebraic variables
     * \a alg and \a z have matching dimensions and
     * <tt>0 == alg(z, ...)</tt> implicitly defines \a z.
     */
    std::vector<MX> z, alg, lam_alg;

    /** \brief Quadrature states
     * Quadrature states are defined by ODEs whose state does not enter in the right-hand-side.
     */
    std::vector<MX> q, quad, lam_quad;


    /** \brief Local variables and corresponding definitions
     */
    std::vector<MX> w, wdef, lam_wdef;

    /** \brief Output variables and corresponding definitions
     */
    std::vector<MX> y, ydef, lam_ydef;

    /** \brief Free controls
     * The trajectories of the free controls are decision variables of the optimal control problem.
     * They are chosen by the optimization algorithm in order to minimize the cost functional.
     */
    std::vector<MX> u;

    /** \brief Parameters
     * A parameter is constant over time, but whose value is chosen by e.g. an
     * optimization algorithm.
     */
    std::vector<MX> p;

    /** \brief Named constants */
    std::vector<MX> c, cdef;

    /** \brief Dependent parameters and corresponding definitions
     * Interdependencies are allowed but must be non-cyclic.
     */
    std::vector<MX> d, ddef, lam_ddef;
    ///@}

    /** \brief Auxiliary variables: Used e.g. to define functions */
    std::vector<MX> aux;

    /** \brief Initial conditions
     * At <tt>t==0</tt>, <tt>0 == init(sdot, s, ...)</tt> holds in addition to
     * the ode and/or dae.
     */
    std::vector<MX> init;
    ///@}

    /** @name Symbolic modeling
     *  Formulate an optimal control problem
     */
    ///@{
    /// Add a new parameter
    MX add_p(const std::string& name=std::string(), casadi_int n=1);

    /// Add a new control
    MX add_u(const std::string& name=std::string(), casadi_int n=1);

    /// Add a new differential state
    MX add_x(const std::string& name=std::string(), casadi_int n=1);

    /// Add a implicit state
    std::pair<MX, MX> add_s(const std::string& name=std::string(), casadi_int n=1);

    /// Add a new algebraic variable
    MX add_z(const std::string& name=std::string(), casadi_int n=1);

    /// Add a new quadrature state
    MX add_q(const std::string& name=std::string(), casadi_int n=1);

    /// Add a new dependent parameter
    MX add_d(const std::string& name, const MX& new_ddef);

    /// Add a new output
    MX add_y(const std::string& name, const MX& new_ydef);

    /// Add an ordinary differential equation
    void add_ode(const std::string& name, const MX& new_ode);

    /// Add a differential-algebraic equation
    void add_dae(const std::string& name, const MX& new_dae);

    /// Add an algebraic equation
    void add_alg(const std::string& name, const MX& new_alg);

    /// Add a quadrature equation
    void add_quad(const std::string& name, const MX& new_quad);

    /// Add an auxiliary variable
    MX add_aux(const std::string& name=std::string(), casadi_int n=1);

    /// Check if dimensions match
    void sanity_check() const;
    ///@}

    /** @name Manipulation
     *  Reformulate the dynamic optimization problem.
     */
    ///@{

    /// Identify and separate the algebraic variables and equations in the DAE
    void split_dae();

    /// Eliminate algebraic variables and equations transforming them into outputs
    void eliminate_alg();

    /// Transform the implicit DAE to a semi-explicit DAE
    void make_semi_explicit();

    /// Transform the implicit DAE or semi-explicit DAE into an explicit ODE
    void make_explicit();

    /// Sort dependent parameters
    void sort_d();

    /// Eliminate interdependencies amongst dependent parameters
    void split_d();

    /// Eliminate dependent parameters
    void eliminate_d();

    /// Eliminate quadrature states and turn them into ODE states
    void eliminate_quad();

    /// Sort the DAE and implicitly defined states
    void sort_dae();

    /// Sort the algebraic equations and algebraic states
    void sort_alg();

    /// Scale the variables
    void scale_variables();

    /// Scale the implicit equations
    void scale_equations();
    ///@}

    /** @name Functions
     *  Add or load auxiliary functions
     */
    ///@{

    /// Add a function from loaded expressions
    Function add_fun(const std::string& name,
                     const std::vector<std::string>& arg,
                     const std::vector<std::string>& res, const Dict& opts=Dict());

    /// Add an already existing function
    Function add_fun(const Function& f);

    /// Add an external function
    Function add_fun(const std::string& name, const Importer& compiler,
                     const Dict& opts=Dict());

    /// Does a particular function already exist?
    bool has_fun(const std::string& name) const;

    /// Get function by name
    Function fun(const std::string& name) const;
  ///@}

    /** @name Import and export
     */
    ///@{
    /// Import existing problem from FMI/XML
    void parse_fmi(const std::string& filename);

#ifndef SWIG
    // Input convension in codegen
    enum DaeBuilderIn {
      DAE_BUILDER_T,
      DAE_BUILDER_C,
      DAE_BUILDER_P,
      DAE_BUILDER_D,
      DAE_BUILDER_U,
      DAE_BUILDER_X,
      DAE_BUILDER_S,
      DAE_BUILDER_SDOT,
      DAE_BUILDER_Z,
      DAE_BUILDER_Q,
      DAE_BUILDER_W,
      DAE_BUILDER_Y,
      DAE_BUILDER_NUM_IN
    };

    // Output convension in codegen
    enum DaeBuilderOut {
      DAE_BUILDER_DDEF,
      DAE_BUILDER_WDEF,
      DAE_BUILDER_ODE,
      DAE_BUILDER_DAE,
      DAE_BUILDER_ALG,
      DAE_BUILDER_QUAD,
      DAE_BUILDER_YDEF,
      DAE_BUILDER_NUM_OUT
    };

    // Get string representation for input, given enum
    static std::string name_in(DaeBuilderIn ind);

    // Get string representation for all inputs
    static std::string name_in();

    // Get enum representation for input, given string
    static DaeBuilderIn enum_in(const std::string& id);

    // Get enum representation for input, given vector of strings
    static std::vector<DaeBuilderIn> enum_in(const std::vector<std::string>& id);

    // Get string representation for output, given enum
    static std::string name_out(DaeBuilderOut ind);

    // Get string representation for all outputs
    static std::string name_out();

    // Get enum representation for output, given string
    static DaeBuilderOut enum_out(const std::string& id);

    // Get enum representation for output, given vector of strings
    static std::vector<DaeBuilderOut> enum_out(const std::vector<std::string>& id);

    // Get input expression, given enum
    std::vector<MX> input(DaeBuilderIn ind) const;

    // Get output expression, given enum
    std::vector<MX> output(DaeBuilderOut ind) const;

    // Get input expression, given enum
    std::vector<MX> input(std::vector<DaeBuilderIn>& ind) const;

    // Get output expression, given enum
    std::vector<MX> output(std::vector<DaeBuilderOut>& ind) const;

    // Get multiplier corresponding to an output expression, given enum
    std::vector<MX> multiplier(DaeBuilderOut ind) const;
#endif // SWIG

    /// Add a named linear combination of output expressions
    MX add_lc(const std::string& name,
              const std::vector<std::string>& f_out);

    /// Construct a function object
    Function create(const std::string& fname,
                    const std::vector<std::string>& s_in,
                    const std::vector<std::string>& s_out) const;
    ///@}

    /// Get variable expression by name
    MX var(const std::string& name) const;

    /// Get variable expression by name
    MX operator()(const std::string& name) const {return var(name);}

    /// Get a derivative expression by name
    MX der(const std::string& name) const;

    /// Get a derivative expression by non-differentiated expression
    MX der(const MX& var) const;

    /// Get the nominal value by name
    double nominal(const std::string& name) const;

    /// Get the nominal value(s) by expression
    std::vector<double> nominal(const MX& var) const;

    /// Set the nominal value by name
    void set_nominal(const std::string& name, double val);

    /// Set the nominal value(s) by expression
    void set_nominal(const MX& var, const std::vector<double>& val);

    /// Get the lower bound by name
    double min(const std::string& name, bool normalized=false) const;

    /// Get the lower bound(s) by expression
    std::vector<double> min(const MX& var, bool normalized=false) const;

    /// Set the lower bound by name
    void set_min(const std::string& name, double val, bool normalized=false);

    /// Set the lower bound(s) by expression
    void set_min(const MX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the upper bound by name
    double max(const std::string& name, bool normalized=false) const;

    /// Get the upper bound(s) by expression
    std::vector<double> max(const MX& var, bool normalized=false) const;

    /// Set the upper bound by name
    void set_max(const std::string& name, double val, bool normalized=false);

    /// Set the upper bound(s) by expression
    void set_max(const MX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the initial guess by name
    double guess(const std::string& name, bool normalized=false) const;

    /// Get the initial guess(es) by expression
    std::vector<double> guess(const MX& var, bool normalized=false) const;

    /// Set the initial guess by name
    void set_guess(const std::string& name, double val, bool normalized=false);

    /// Set the initial guess(es) by expression
    void set_guess(const MX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) value at time 0 by name
    double start(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) value(s) at time 0 by expression
    std::vector<double> start(const MX& var, bool normalized=false) const;

    /// Set the (optionally normalized) value at time 0 by name
    void set_start(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) value(s) at time 0 by expression
    void set_start(const MX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) derivative value at time 0 by name
    double derivative_start(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) derivative value(s) at time 0 by expression
    std::vector<double> derivative_start(const MX& var, bool normalized=false) const;

    /// Set the (optionally normalized) derivative value at time 0 by name
    void set_derivative_start(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) derivative value(s) at time 0 by expression
    void set_derivative_start(const MX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the unit for a component
    std::string unit(const std::string& name) const;

    /// Get the unit given a vector of symbolic variables (all units must be identical)
    std::string unit(const MX& var) const;

    /// Set the unit for a component
    void set_unit(const std::string& name, const std::string& val);

    /// Readable name of the class
    std::string type_name() const {return "DaeBuilder";}

    ///  Print representation
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const {
      std::stringstream ss;
      disp(ss, more);
      return ss.str();
    }

    /// Add a variable
    void add_variable(const std::string& name, const Variable& var);

    /// Add a new variable: returns corresponding symbolic expression
    MX add_variable(const std::string& name, casadi_int n=1);

    /// Add a new variable: returns corresponding symbolic expression
    MX add_variable(const std::string& name, const Sparsity& sp);

    ///@{
    /// Access a variable by name
    Variable& variable(const std::string& name);
    const Variable& variable(const std::string& name) const;
    ///@}

#ifndef SWIG
    // Internal methods
  protected:

    /// Get the qualified name
    static std::string qualified_name(const XmlNode& nn);

    /// Find of variable by name
    typedef std::map<std::string, Variable> VarMap;
    VarMap varmap_;

    /// Linear combinations of output expressions
    std::map<std::string, MX> lin_comb_;

    /** \brief Functions */
    std::vector<Function> fun_;

    /// Read an equation
    MX read_expr(const XmlNode& odenode);

    /// Read a variable
    Variable& read_variable(const XmlNode& node);

    /// Get an attribute by expression
    typedef double (DaeBuilder::*getAtt)(const std::string& name, bool normalized) const;
    std::vector<double> attribute(getAtt f, const MX& var, bool normalized) const;

    /// Get a symbolic attribute by expression
    typedef MX (DaeBuilder::*getAttS)(const std::string& name) const;
    MX attribute(getAttS f, const MX& var) const;

    /// Set an attribute by expression
    typedef void (DaeBuilder::*setAtt)(const std::string& name, double val, bool normalized);
    void set_attribute(setAtt f, const MX& var, const std::vector<double>& val, bool normalized);

    /// Set a symbolic attribute by expression
    typedef void (DaeBuilder::*setAttS)(const std::string& name, const MX& val);
    void set_attribute(setAttS f, const MX& var, const MX& val);

#endif // SWIG

  };

} // namespace casadi

#endif // CASADI_DAE_BUILDER_HPP
