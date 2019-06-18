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

#ifndef CASADI_OPTISTACK_INTERNAL_HPP
#define CASADI_OPTISTACK_INTERNAL_HPP

#include "optistack.hpp"
#include "shared_object_internal.hpp"

namespace casadi {

#ifndef SWIG
/// Pointer that gets set to null when copied
template<class T>
class null_ptr_on_copy {
public:
  null_ptr_on_copy<T>() : ptr_(nullptr) {}
  null_ptr_on_copy<T>(const null_ptr_on_copy<T>& rhs) : ptr_(nullptr) {}
  void operator=(T* ptr) { ptr_ = ptr; }
  T* operator->() { return ptr_; }
  operator bool() const { return ptr_; }
private:
  T* ptr_;
};
#endif

/** \brief A simplified interface for NLP modeling/solving


      \date 2017
      \author Joris Gillis, Erik Lambrechts
*/
class CASADI_EXPORT OptiNode :
    public SharedObjectInternal {
  friend class InternalOptiCallback;
public:

  /// Create Opti Context
  OptiNode();

  /// Destructor
  ~OptiNode();

  /// Create a decision variable (symbol)
  MX variable(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");

  /// Create a parameter (symbol); fixed during optimization
  MX parameter(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");

  /// Set objective
  void minimize(const MX& f);

  /// brief Add constraints
  void subject_to(const MX& g);
  /// Clear constraints
  void subject_to();

  /// Solver
  void solver(const std::string& solver,
              const Dict& plugin_options=Dict(),
              const Dict& solver_options=Dict());

  /// @{
  /// Set initial value for decision variables
  void set_initial(const MX& x, const DM& v);
  void set_initial(const std::vector<MX>& assignments);
  /// @}

  /// @{
  /** \brief Set value of parameter
  *
  * Each parameter must be given a value before 'solve' can be called
  */
  void set_value(const MX& x, const DM& v);
  void set_value(const std::vector<MX>& assignments);
  /// @}

  /// Crunch the numbers; solve the problem
  OptiSol solve();

  /// @{
  /// Obtain value of expression at the current value
  DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
  DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const { return x; }
  DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const {
    return DM::nan(x.sparsity());
  }
  /// @}

  /// Copy
  Opti copy() const;

  /// Get statistics
  Dict stats() const;

  /// Get return status of solver
  std::string return_status() const;
  /// Did the solver return successfully?
  bool return_success() const;

  /// Get the underlying CasADi solver of the Opti stack
  Function casadi_solver() const;

  /// get assignment expressions for initial values
  std::vector<MX> initial() const;

  /// get assignment expressions for latest values
  std::vector<MX> value_variables() const;
  std::vector<MX> value_parameters() const;

  void callback_class(OptiCallback* callback);
  void callback_class();
  bool has_callback_class() const;

  /// return true if expression is only dependant on Opti parameters, not variables
  bool is_parametric(const MX& expr) const;

  /// @{
  /// Get symbols present in expression
  std::vector<MX> symvar() const;
  std::vector<MX> symvar(const MX& expr) const;
  std::vector<MX> symvar(const MX& expr, VariableType type) const;
  /// @}

  /// Interpret an expression (for internal use only)
  MetaCon canon_expr(const MX& expr) const;

  /// Get meta-data of symbol (for internal use only)
  MetaVar get_meta(const MX& m) const;

  /// Get meta-data of symbol (for internal use only)
  MetaCon get_meta_con(const MX& m) const;

  /// Set meta-data of an expression
  void set_meta(const MX& m, const MetaVar& meta);

  /// Set meta-data of an expression
  void set_meta_con(const MX& m, const MetaCon& meta);

  /// add meta-data of an expression
  void update_user_dict(const MX& m, const Dict& meta);
  Dict user_dict(const MX& m) const;

  /// get the dual variable
  MX dual(const MX& m) const;

  void assert_active_symbol(const MX& m) const;

  std::vector<MX> active_symvar(VariableType type) const;
  std::vector<DM> active_values(VariableType type) const;

  MX x_lookup(casadi_int i) const;
  MX g_lookup(casadi_int i) const;

  std::string x_describe(casadi_int i) const;
  std::string g_describe(casadi_int i) const;
  std::string describe(const MX& x, casadi_int indent=0) const;

  void solve_prepare();
  DMDict solve_actual(const DMDict& args);

  DMDict arg() const { return arg_; }
  void res(const DMDict& res);
  DMDict res() const { return res_; }
  std::vector<MX> constraints() const { return g_; }
  MX objective() const { return f_; }

  OptiAdvanced baked_copy() const {
    OptiAdvanced s = copy();
    if (s.problem_dirty()) s.bake();
    return s;
  }

  std::string class_name() const override { return "OptiNode"; }

  /// Number of (scalarised) decision variables
  casadi_int nx() const {
    if (problem_dirty()) return baked_copy().nx();
    return nlp_.at("x").size1();
  }

  /// Number of (scalarised) parameters
  casadi_int np() const {
    if (problem_dirty()) return baked_copy().np();
    return nlp_.at("p").size1();
  }

  /// Number of (scalarised) constraints
  casadi_int ng() const {
    if (problem_dirty()) return baked_copy().ng();
    return nlp_.at("g").size1();
  }

  /// Get all (scalarised) decision variables as a symbolic column vector
  MX x() const {
    if (problem_dirty()) return baked_copy().x();
    return nlp_.at("x");
  }

  /// Get all (scalarised) parameters as a symbolic column vector
  MX p() const {
    if (problem_dirty()) return baked_copy().p();
    return nlp_.at("p");
  }

  /// Get all (scalarised) constraint expressions as a column vector
  MX g() const {
    if (problem_dirty()) return baked_copy().g();
    return nlp_.at("g");
  }

  /// Get objective expression
  MX f() const {
    if (problem_dirty()) return baked_copy().f();
    return nlp_.at("f");
  }

  MX lbg() const {
    if (problem_dirty()) return baked_copy().lbg();
    return bounds_lbg_;
  }

  MX ubg() const {
    if (problem_dirty()) return baked_copy().ubg();
    return bounds_ubg_;
  }

  /// Get dual variables as a symbolic column vector
  MX lam_g() const {
    if (problem_dirty()) return baked_copy().lam_g();
    return lam_;
  }
  void assert_empty() const;

  ///  Print representation
  void disp(std::ostream& stream, bool more=false) const override;

  /// Fix the structure of the optimization problem
  void bake();

  casadi_int instance_number() const;

  static OptiNode* create();

  bool problem_dirty_;
  void mark_problem_dirty(bool flag=true) { problem_dirty_=flag; mark_solver_dirty(); }
  bool problem_dirty() const { return problem_dirty_; }

  bool solver_dirty_;
  void mark_solver_dirty(bool flag=true) { solver_dirty_=flag; mark_solved(false); }
  bool solver_dirty() const { return solver_dirty_; }

  bool solved_;
  void mark_solved(bool flag=true) { solved_ = flag;}
  bool solved() const { return solved_; }

  void assert_solved() const;
  void assert_baked() const;

private:

  static std::map<VariableType, std::string> VariableType2String_;
  std::string variable_type_to_string(VariableType vt) const;

  bool parse_opti_name(const std::string& name, VariableType& vt) const;
  void register_dual(MetaCon& meta);

  /// Set value of symbol
  void set_value_internal(const MX& x, const DM& v);

  /** \brief decompose a chain of inequalities
  *
  * a<=b -> [a,b]
  * a<=b<=c [a,b,c]
  *
  * When flipped is set, [a,b,c] corresponds to c>=b>=a
  */
  static std::vector<MX> ineq_unchain(const MX& a, bool& SWIG_OUTPUT(flipped));

  /// Get meta-dat by const-ref
  const MetaVar& meta(const MX& m) const;
  /// Get meta-dat by ref
  MetaVar& meta(const MX& m);

  /// Get meta-dat by const-ref
  const MetaCon& meta_con(const MX& m) const;
  /// Get meta-dat by ref
  MetaCon& meta_con(const MX& m);

  /// Sort symbols according to Opti order
  std::vector<MX> sort(const std::vector<MX>& v) const;

  /// Throw an error
  void assert_has(const MX& m) const;

  bool has(const MX& m) const;

  /// Throw an error
  void assert_has_con(const MX& m) const;

  bool has_con(const MX& m) const;

  // Data members

  /// Map symbols to metadata
  std::map<MXNode*, MetaVar> meta_;
  /// map constraints to metadata
  std::map<MXNode*, MetaCon> meta_con_;

  /// Store references to all symbols
  std::vector<MX> symbols_;

  /// Symbol counter
  casadi_int count_;

  casadi_int count_var_;
  casadi_int count_par_;
  casadi_int count_dual_;

  /// Storing initial/latest values for all variables (including inactive)
  std::map< VariableType, std::vector<DM> > store_initial_, store_latest_;

  /// Is symbol present in problem?
  std::vector<bool> symbol_active_;

  /// Solver
  Function solver_;

  /// Result of solver
  DMDict res_;
  DMDict arg_;
  MXDict nlp_;
  MX lam_;

  /// Bounds helper function: p -> lbg, ubg
  Function bounds_;
  MX bounds_lbg_;
  MX bounds_ubg_;

  /// Constraints verbatim as passed in with 'subject_to'
  std::vector<MX> g_;

  /// Objective verbatim as passed in with 'minimize'
  MX f_;

  null_ptr_on_copy<OptiCallback> user_callback_;
  Function callback_;

  bool old_callback() const;

  std::string solver_name_;
  Dict solver_options_;

  void assert_only_opti_symbols(const MX& e) const;
  void assert_only_opti_nondual(const MX& e) const;


  static casadi_int instance_count_;
  casadi_int instance_number_;


  std::string name_prefix() const;

  static std::string format_stacktrace(const Dict& stacktrace, casadi_int indent);

};


} // namespace casadi

#endif // CASADI_OPTISTACK_INTERNAL_HPP
