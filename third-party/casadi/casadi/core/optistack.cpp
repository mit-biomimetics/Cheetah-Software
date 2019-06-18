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

#include "optistack_internal.hpp"
#include "exception.hpp"
#include "global_options.hpp"

using namespace std;
namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in Opti::" FNAME " "\
  "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
  + string(WHAT));


OptiNode* Opti::operator->() {
  if (getCount()!=1) *this = copy();
  return static_cast<OptiNode*>(SharedObject::operator->());
}

const OptiNode* Opti::operator->() const {
  return static_cast<const OptiNode*>(SharedObject::operator->());
}


Opti::Opti() {
  own(OptiNode::create());
}

MX Opti::variable(casadi_int n, casadi_int m, const std::string& attribute) {
  try {
    return (*this)->variable(n, m, attribute);
  } catch (exception& e) {
    THROW_ERROR("variable", e.what());
  }
}

Opti::Opti(OptiNode* node) {
  own(node);
}

Opti::Opti(const Opti& rhs) : SharedObject(rhs) {
  callback_class();
}

OptiAdvanced::OptiAdvanced(const Opti& rhs) : Opti(rhs) {
}

Opti Opti::create(OptiNode* node) {
  return Opti(node);
}

MX Opti::parameter(casadi_int n, casadi_int m, const std::string& attribute) {
  try {
    return (*this)->parameter(n, m, attribute);
  } catch (exception& e) {
    THROW_ERROR("parameter", e.what());
  }
}

void Opti::minimize(const MX& f) {
  try {
    (*this)->minimize(f);
  } catch (exception& e) {
    THROW_ERROR("minimize", e.what());
  }
}

void Opti::subject_to(const MX& g) {
  try {
    (*this)->subject_to(g);
  } catch (exception& e) {
    THROW_ERROR("subject_to", e.what());
  }
}

void Opti::subject_to(const std::vector<MX>& g) {
  for (const auto& gs : g) subject_to(gs);
}

void Opti::subject_to() {
  try {
    (*this)->subject_to();
  } catch (exception& e) {
    THROW_ERROR("subject_to", e.what());
  }
}


void Opti::solver(const std::string& solver,
                       const Dict& plugin_options,
                       const Dict& solver_options) {
  try {
    (*this)->solver(solver, plugin_options, solver_options);
  } catch (exception& e) {
    THROW_ERROR("solver", e.what());
  }
}

void Opti::set_initial(const MX& x, const DM& v) {
  try {
    (*this)->set_initial(x, v);
  } catch (exception& e) {
    THROW_ERROR("set_initial", e.what());
  }
}
void Opti::set_initial(const std::vector<MX>& assignments) {
  try {
    (*this)->set_initial(assignments);
  } catch (exception& e) {
    THROW_ERROR("set_initial", e.what());
  }
}


void Opti::set_value(const MX& x, const DM& v) {
  try {
    (*this)->set_value(x, v);
  } catch (exception& e) {
    THROW_ERROR("set_value", e.what());
  }
}

void Opti::set_value(const std::vector<MX>& assignments) {
  try {
    (*this)->set_value(assignments);
  } catch (exception& e) {
    THROW_ERROR("set_value", e.what());
  }
}

OptiSol Opti::solve() {
  try {
    return (*this)->solve();
  } catch (exception& e) {
    THROW_ERROR("solve", e.what());
  }
}

DM Opti::value(const MX& x, const std::vector<MX>& values) const {
  try {
    return (*this)->value(x, values);
  } catch (exception& e) {
    THROW_ERROR("value", e.what());
  }
}


DM Opti::value(const DM& x, const std::vector<MX>& values) const {
  try {
    return (*this)->value(x, values);
  } catch (exception& e) {
    THROW_ERROR("value", e.what());
  }
}

DM Opti::value(const SX& x, const std::vector<MX>& values) const {
  try {
    return (*this)->value(x, values);
  } catch (exception& e) {
    THROW_ERROR("value", e.what());
  }
}

Dict Opti::stats() const {
  try {
    return (*this)->stats();
  } catch (exception& e) {
    THROW_ERROR("stats", e.what());
  }
}

std::string Opti::return_status() const {
  try {
    return (*this)->return_status();
  } catch (exception& e) {
    THROW_ERROR("return_status", e.what());
  }
}

std::vector<MX> Opti::initial() const {
  try {
    return (*this)->initial();
  } catch (exception& e) {
    THROW_ERROR("initial", e.what());
  }
}

std::vector<MX> Opti::value_variables() const {
  try {
    return (*this)->value_variables();
  } catch (exception& e) {
    THROW_ERROR("value_variables", e.what());
  }
}

std::vector<MX> Opti::value_parameters() const {
  try {
    return (*this)->value_parameters();
  } catch (exception& e) {
    THROW_ERROR("value_parameters", e.what());
  }
}

MX Opti::dual(const MX& m) const {
  try {
    return (*this)->dual(m);
  } catch (exception& e) {
    THROW_ERROR("dual", e.what());
  }
}

casadi_int Opti::nx() const {
  try {
    return (*this)->nx();
  } catch (exception& e) {
    THROW_ERROR("nx", e.what());
  }
}

casadi_int Opti::np() const {
  try {
    return (*this)->np();
  } catch (exception& e) {
    THROW_ERROR("nx", e.what());
  }
}

casadi_int Opti::ng() const {
  try {
    return (*this)->ng();
  } catch (exception& e) {
    THROW_ERROR("ng", e.what());
  }
}

MX Opti::x() const {
  try {
    return (*this)->x();
  } catch (exception& e) {
    THROW_ERROR("x", e.what());
  }
}

MX Opti::p() const {
  try {
    return (*this)->p();
  } catch (exception& e) {
    THROW_ERROR("p", e.what());
  }
}

MX Opti::g() const {
  try {
    return (*this)->g();
  } catch (exception& e) {
    THROW_ERROR("g", e.what());
  }
}

MX Opti::f() const {
  try {
    return (*this)->f();
  } catch (exception& e) {
    THROW_ERROR("f", e.what());
  }
}

MX Opti::lbg() const {
  try {
    return (*this)->lbg();
  } catch (exception& e) {
    THROW_ERROR("lbg", e.what());
  }
}

MX Opti::ubg() const {
  try {
    return (*this)->ubg();
  } catch (exception& e) {
    THROW_ERROR("ubg", e.what());
  }
}


MX Opti::lam_g() const {
  try {
    return (*this)->lam_g();
  } catch (exception& e) {
    THROW_ERROR("lam_g", e.what());
  }
}

void Opti::callback_class(OptiCallback* callback) {
  try {
    (*this)->callback_class(callback);
  } catch (exception& e) {
    THROW_ERROR("callback_class", e.what());
  }
}

void Opti::callback_class() {
  try {
    if ((*this)->has_callback_class() && getCount()!=1) {
      Opti ret = copy();
      ret.callback_class();
      *this = ret;
    } else {
      (*this)->callback_class();
    }
  } catch (exception& e) {
    THROW_ERROR("callback_class", e.what());
  }
}

void Opti::update_user_dict(const MX& m, const Dict& meta) {
  try {
    (*this)->update_user_dict(m, meta);
  } catch (exception& e) {
    THROW_ERROR("update_user_dict", e.what());
  }
}

void Opti::update_user_dict(const std::vector<MX>& ms, const Dict& meta) {
  for (const auto& m : ms)
     update_user_dict(m, meta);
}

Dict Opti::user_dict(const MX& m) const {
  try {
    return (*this)->user_dict(m);
  } catch (exception& e) {
    THROW_ERROR("user_dict", e.what());
  }
}

Function OptiAdvanced::casadi_solver() const {
  try {
    return (*this)->casadi_solver();
  } catch (exception& e) {
    THROW_ERROR("casadi_solver", e.what());
  }
}

bool OptiAdvanced::is_parametric(const MX& expr) const {
  try {
    return (*this)->is_parametric(expr);
  } catch (exception& e) {
    THROW_ERROR("is_parametric", e.what());
  }
}

std::vector<MX> OptiAdvanced::symvar() const {
  try {
    return (*this)->symvar();
  } catch (exception& e) {
    THROW_ERROR("symvar", e.what());
  }
}

std::vector<MX> OptiAdvanced::symvar(const MX& expr) const {
  try {
    return (*this)->symvar(expr);
  } catch (exception& e) {
    THROW_ERROR("symvar", e.what());
  }
}

std::vector<MX> OptiAdvanced::symvar(const MX& expr, VariableType type) const {
  try {
    return (*this)->symvar(expr, type);
  } catch (exception& e) {
    THROW_ERROR("symvar", e.what());
  }
}

MetaCon OptiAdvanced::canon_expr(const MX& expr) const {
  try {
    return (*this)->canon_expr(expr);
  } catch (exception& e) {
    THROW_ERROR("canon_expr", e.what());
  }
}

MetaVar OptiAdvanced::get_meta(const MX& m) const {
  try {
    return (*this)->get_meta(m);
  } catch (exception& e) {
    THROW_ERROR("get_meta", e.what());
  }
}

MetaCon OptiAdvanced::get_meta_con(const MX& m) const {
  try {
    return (*this)->get_meta_con(m);
  } catch (exception& e) {
    THROW_ERROR("get_meta_con", e.what());
  }
}

void OptiAdvanced::set_meta(const MX& m, const MetaVar& meta) {
  try {
    (*this)->set_meta(m, meta);
  } catch (exception& e) {
    THROW_ERROR("set_meta", e.what());
  }
}

void OptiAdvanced::set_meta_con(const MX& m, const MetaCon& meta) {
  try {
    return (*this)->set_meta_con(m, meta);
  } catch (exception& e) {
    THROW_ERROR("set_meta_con", e.what());
  }
}


void OptiAdvanced::assert_active_symbol(const MX& m) const {
  try {
    (*this)->assert_active_symbol(m);
  } catch (exception& e) {
    THROW_ERROR("assert_active_symbol", e.what());
  }
}

std::vector<MX> OptiAdvanced::active_symvar(VariableType type) const {
  try {
    return (*this)->active_symvar(type);
  } catch (exception& e) {
    THROW_ERROR("active_symvar", e.what());
  }
}

std::vector<DM> OptiAdvanced::active_values(VariableType type) const {
  try {
    return (*this)->active_values(type);
  } catch (exception& e) {
    THROW_ERROR("active_values", e.what());
  }
}

MX OptiAdvanced::x_lookup(casadi_int i) const {
  try {
    return (*this)->x_lookup(i);
  } catch (exception& e) {
    THROW_ERROR("x_lookup", e.what());
  }
}

MX OptiAdvanced::g_lookup(casadi_int i) const {
  try {
    return (*this)->g_lookup(i);
  } catch (exception& e) {
    THROW_ERROR("g_lookup", e.what());
  }
}

std::string OptiAdvanced::x_describe(casadi_int i) const {
  try {
    return (*this)->x_describe(i);
  } catch (exception& e) {
    THROW_ERROR("x_describe", e.what());
  }
}
std::string OptiAdvanced::g_describe(casadi_int i) const {
  try {
    return (*this)->g_describe(i);
  } catch (exception& e) {
    THROW_ERROR("g_describe", e.what());
  }
}
std::string OptiAdvanced::describe(const MX& x, casadi_int indent) const {
  try {
    return (*this)->describe(x, indent);
  } catch (exception& e) {
    THROW_ERROR("describe", e.what());
  }
}

void OptiAdvanced::show_infeasibilities(double tol) const {
  std::vector<double> g_ = value(g()).get_elements();
  std::vector<double> lbg_ = value(lbg()).get_elements();
  std::vector<double> ubg_ = value(ubg()).get_elements();

  uout() << "Violated constraints (tol " << tol << "), in order of declaration:" << std::endl;
  for (casadi_int i=0;i<g_.size();++i) {
    double err = std::max(g_[i]-ubg_[i], lbg_[i]-g_[i]);
    if (err>=tol) {
      uout() << "------- i = " << i+GlobalOptions::start_index;
      uout() << "/" << g_.size() << " ------ " << std::endl;
      uout() << lbg_[i] << " <= " << g_[i] << " <= " << ubg_[i];
      uout() << " (viol " << err << ")" << std::endl;
      uout() << g_describe(i) << std::endl;
    }
  }
}

void OptiAdvanced::solve_prepare() {
  try {
    (*this)->solve_prepare();
  } catch (exception& e) {
    THROW_ERROR("solve_prepare", e.what());
  }
}
DMDict OptiAdvanced::solve_actual(const DMDict& args) {
  try {
    return (*this)->solve_actual(args);
  } catch (exception& e) {
    THROW_ERROR("solve_actual", e.what());
  }
}

DMDict OptiAdvanced::arg() const {
  try {
    return (*this)->arg();
  } catch (exception& e) {
    THROW_ERROR("arg", e.what());
  }
}


void OptiAdvanced::res(const DMDict& res) {
  try {
    return (*this)->res(res);
  } catch (exception& e) {
    THROW_ERROR("res", e.what());
  }
}

DMDict OptiAdvanced::res() const {
  try {
    return (*this)->res();
  } catch (exception& e) {
    THROW_ERROR("res", e.what());
  }
}

std::vector<MX> OptiAdvanced::constraints() const {
  try {
    return (*this)->constraints();
  } catch (exception& e) {
    THROW_ERROR("constraints", e.what());
  }
}
MX OptiAdvanced::objective() const {
  return f();
}

OptiAdvanced OptiAdvanced::baked_copy() const {
  try {
    return (*this)->baked_copy();
  } catch (exception& e) {
    THROW_ERROR("baked_copy", e.what());
  }
}

void OptiAdvanced::assert_empty() const {
  try {
    return (*this)->assert_empty();
  } catch (exception& e) {
    THROW_ERROR("assert_empty", e.what());
  }
}

casadi_int OptiAdvanced::instance_number() const {
  try {
    return (*this)->instance_number();
  } catch (exception& e) {
    THROW_ERROR("instance_number", e.what());
  }
}

void Opti::disp(std::ostream& stream, bool more) const {
  stream << "Opti {" << std::endl;
  OptiAdvanced mycopy = debug();
  stream << "  instance #" << mycopy.instance_number() << std::endl;
  if (mycopy.problem_dirty()) mycopy.bake();
  stream << "  #variables: " << mycopy.active_symvar(OPTI_VAR).size()
    << " (nx = " << mycopy.nx() << ")" <<  std::endl;
  stream << "  #parameters: " << mycopy.active_symvar(OPTI_PAR).size()
    << " (np = " << mycopy.np() << ")" << std::endl;
  stream << "  #constraints: " << mycopy.active_symvar(OPTI_DUAL_G).size()
    << " (ng = " << mycopy.ng() << ")" << std::endl;
  if (mycopy.solver_dirty()) {
    stream << "  CasADi solver needs updating." << std::endl;
  } else {
    stream << "  CasADi solver allocated." << std::endl;
  }
  if (mycopy.solved()) {
    stream << "  CasADi solver was called: " << mycopy.return_status() << std::endl;
  }
  stream << "}";
}

std::string Opti::get_str(bool more) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
}

void OptiAdvanced::bake() {
  try {
    (*this)->bake();
  } catch (exception& e) {
    THROW_ERROR("bake", e.what());
  }
}

void OptiAdvanced::mark_problem_dirty(bool flag) {
  try {
    (*this)->mark_problem_dirty(flag);
  } catch (exception& e) {
    THROW_ERROR("mark_problem_dirty", e.what());
  }
}
bool OptiAdvanced::problem_dirty() const {
  try {
    return (*this)->problem_dirty();
  } catch (exception& e) {
    THROW_ERROR("problem_dirty", e.what());
  }
}

void OptiAdvanced::mark_solver_dirty(bool flag) {
  try {
    (*this)->mark_solver_dirty(flag);
  } catch (exception& e) {
    THROW_ERROR("mark_solver_dirty", e.what());
  }
}
bool OptiAdvanced::solver_dirty() const {
  try {
    return (*this)->solver_dirty();
  } catch (exception& e) {
    THROW_ERROR("solver_dirty", e.what());
  }
}

void OptiAdvanced::mark_solved(bool flag) {
  try {
    (*this)->mark_solved(flag);
  } catch (exception& e) {
    THROW_ERROR("mark_solved", e.what());
  }
}

bool OptiAdvanced::solved() const {
  try {
    return (*this)->solved();
  } catch (exception& e) {
    THROW_ERROR("solved", e.what());
  }
}

void OptiAdvanced::assert_solved() const {
  try {
    (*this)->assert_solved();
  } catch (exception& e) {
    THROW_ERROR("assert_solved", e.what());
  }
}
void OptiAdvanced::assert_baked() const {
  try {
    (*this)->assert_baked();
  } catch (exception& e) {
    THROW_ERROR("assert_baked", e.what());
  }
}

OptiAdvanced Opti::debug() const {
  return *this;
}
OptiAdvanced Opti::advanced() const {
  return *this;
}
Opti Opti::copy() const {
  return (*this)->copy();
}

OptiSol::OptiSol(const Opti& opti) : optistack_(opti) {
}

void OptiSol::disp(std::ostream& stream, bool more) const {
  optistack_.disp(stream, more);
}

std::string OptiSol::get_str(bool more) const {
  return optistack_.get_str(more);
}

DM OptiSol::value(const MX& x, const std::vector<MX>& values) const {
  return optistack_.value(x, values);
}
DM OptiSol::value(const DM& x, const std::vector<MX>& values) const {
  return optistack_.value(x, values);
}
DM OptiSol::value(const SX& x, const std::vector<MX>& values) const {
  return optistack_.value(x, values);
}

std::vector<MX> OptiSol::value_variables() const {
  return optistack_.value_variables();
}

std::vector<MX> OptiSol::value_parameters() const {
  return optistack_.value_parameters();
}

Dict OptiSol::stats() const {
  return optistack_.stats();
}




} // namespace casadi
