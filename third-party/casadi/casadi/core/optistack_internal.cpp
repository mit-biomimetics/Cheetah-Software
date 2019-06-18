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
#include "nlpsol.hpp"
#include "function_internal.hpp"
#include "global_options.hpp"

using namespace std;
namespace casadi {

class InternalOptiCallback : public FunctionInternal {
  public:

  InternalOptiCallback(OptiNode& sol) : FunctionInternal(class_name()), sol_(sol) {}

  ~InternalOptiCallback() override {}

  /** \brief Get type name */
  std::string class_name() const override {return "InternalOptiCallback";}

  // Number of inputs and outputs
  size_t get_n_in() override { return nlpsol_n_out();}

  Sparsity get_sparsity_in(casadi_int i) override {
    std::string n = nlpsol_out(i);
    casadi_int size = 0;
    if (n=="f") {
      size = 1;
    } else if (n=="lam_x" || n=="x") {
      size = sol_.nx();
    } else if (n=="lam_g" || n=="g") {
      size = sol_.ng();
    } else if (n=="p" || n=="lam_p") {
      size = sol_.np();
      if (size==0) return Sparsity::dense(0, 0);
    } else {
      return Sparsity::dense(0, 0);
    }
    return Sparsity::dense(size, 1);
  }

  void reset() { i=0; }

  /// Evaluate the function numerically
  std::vector<DM> eval_dm(const std::vector<DM>& arg) const override {
    DMDict r;

    for (casadi_int i=0;i<nlpsol_n_out();++i) {
      r[nlpsol_out(i)] = arg[i];
    }

    sol_.res(r);

    if (sol_.user_callback_) sol_.user_callback_->call(i);

    i+=1;
    return {0};
  }

  bool associated_with(const OptiNode* o) { return &sol_==o; }

  private:
    OptiNode& sol_;
    mutable casadi_int i;
};

OptiNode* OptiNode::create() {
return new OptiNode();
}


void OptiNode::callback_class(OptiCallback* callback) {
  user_callback_ = callback;
}

void OptiNode::callback_class() {
  user_callback_ = nullptr;
}

bool OptiNode::has_callback_class() const {
  return user_callback_ != 0;
}

std::string OptiNode::format_stacktrace(const Dict& stacktrace, casadi_int indent) {
  std::string s_indent;
  for (casadi_int i=0;i<indent;++i) {
    s_indent+= "  ";
  }
  std::string description;
  std::string filename = stacktrace.at("file").as_string();
  casadi_int line = stacktrace.at("line").as_int();
  description += "defined at " + filename +":"+str(line);
  std::string name = stacktrace.at("name").as_string();
  if (name!="Unknown" && name!= "<module>")
    description += " in " + stacktrace.at("name").as_string();
  try {
    ifstream file(filename);
    for (casadi_int i=0; i<line-1; ++i) {
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::string contents; std::getline(file, contents);
    auto it = contents.find_first_not_of(" \n");
    if (it!=std::string::npos) {
      description += "\n" + s_indent + contents.substr(it);
    }
  } catch(...) {
    // pass
  }
  return description;
}

std::string OptiNode::describe(const MX& expr, casadi_int indent) const {
  if (problem_dirty()) return baked_copy().describe(expr, indent);
  std::string s_indent;
  for (casadi_int i=0;i<indent;++i) {
    s_indent+= "  ";
  }
  std::string description = s_indent;
  if (expr.is_symbolic()) {
    if (has(expr)) {
      description += "Opti " + variable_type_to_string(meta(expr).type) + " '" + expr.name() +
        "' of shape " + expr.dim();
      const Dict& extra = meta(expr).extra;
      auto it = extra.find("stacktrace");
      if (it!=extra.end()) {
        const Dict& stacktrace = it->second.as_dict();
        description += ", " + format_stacktrace(stacktrace, indent+1);
      }
    } else {
      VariableType vt;
      if (parse_opti_name(expr.name(), vt)) {
        description += "Opti " + variable_type_to_string(vt) + " '" + expr.name() +
          "' of shape " + expr.dim()+
          ", belonging to a different instance of Opti.";
      } else {
        description += "MX symbol '" + expr.name() + "' of shape " + expr.dim();
        description += ", declared outside of Opti.";
      }
    }
  } else {
    if (has_con(expr)) {
      description = "Opti constraint of shape " + expr.dim();
      const Dict& extra = meta_con(expr).extra;
      auto it = extra.find("stacktrace");
      if (it!=extra.end()) {
        const Dict& stacktrace = it->second.as_dict();
        description += ", " + format_stacktrace(stacktrace, indent+1);
      }
    } else {
      std::vector<MX> s = symvar(expr);
      if (s.size()==0) {
        description+= "Constant epxression.";
      } else {
        description+= "General expression, dependent on " + str(s.size()) + " symbols:";
        for (casadi_int i=0;i<s.size();++i) {
          description+= "\n"+describe(s[i], indent+1);
          if (i>5) {
            description+= "\n...";
            break;
          }
        }
      }
    }
  }

  return description;
}

std::string OptiNode::g_describe(casadi_int i) const {
  if (problem_dirty()) return baked_copy().g_describe(i);
  MX expr = g_lookup(i);
  casadi_int local_i = i-meta_con(expr).start + GlobalOptions::start_index;
  std::string description = describe(expr);
  if (expr.numel()>1)
    description += "\nAt nonzero " + str(local_i) + ".";
  return description;
}

std::string OptiNode::x_describe(casadi_int i) const {
  if (problem_dirty()) return baked_copy().x_describe(i);
  MX symbol = x_lookup(i);
  casadi_int local_i = i-meta(symbol).start + GlobalOptions::start_index;
  std::string description = describe(symbol);
  if (symbol.numel()>1)
    description += "\nAt nonzero " + str(local_i) + ".";
  return description;
}

MX OptiNode::x_lookup(casadi_int i) const {
  if (problem_dirty()) return baked_copy().x_lookup(i);
  casadi_assert_dev(i>=0);
  casadi_assert_dev(i<nx());
  std::vector<MX> x = active_symvar(OPTI_VAR);
  for (const auto& e : x) {
    const MetaVar& m = meta(e);
    if (i>=m.start && i<m.stop) return e;
  }
  casadi_error("Internal error");
  return MX();
}

MX OptiNode::g_lookup(casadi_int i) const {
  if (problem_dirty()) return baked_copy().g_lookup(i);
  casadi_assert_dev(i>=0);
  casadi_assert_dev(i<ng());
  for (const auto& e : g_) {
    const MetaCon& m = meta_con(e);
    if (i>=m.start && i<m.stop) return e;
  }
  casadi_error("Internal error");
  return MX();
}

OptiNode::OptiNode() : count_(0), count_var_(0), count_par_(0), count_dual_(0) {
  f_ = 0;
  instance_number_ = instance_count_++;
  user_callback_ = nullptr;
  store_initial_[OPTI_VAR] = {};
  store_initial_[OPTI_PAR] = {};
  store_initial_[OPTI_DUAL_G] = {};
  store_latest_[OPTI_VAR] = {};
  store_latest_[OPTI_DUAL_G] = {};
  mark_problem_dirty();
}

OptiNode::~OptiNode() {
}

MX OptiNode::variable(casadi_int n, casadi_int m, const std::string& attribute) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_VAR;
  meta_data.count = count_++;
  meta_data.i = count_var_++;

  MX symbol, ret;

  if (attribute=="symmetric") {
    casadi_assert(n==m, "You specified attribute 'symmetric', "
      "while matrix is not even square, but " + str(n) + "-by-" + str(m) + ".");
    symbol = MX::sym(name_prefix() + "x_" + str(count_var_), n*(n+1)/2);
    ret = tril2symm(MX(Sparsity::lower(n), symbol));
  } else if (attribute=="full") {
    symbol = MX::sym(name_prefix() + "x_" + str(count_var_), n, m);
    ret = symbol;
  } else {
    casadi_error("Unknown attribute '" + attribute + "'. Choose from 'full' or 'symmetric'.");
  }

  // Store the symbol; preventing it from going ut of scope
  symbols_.push_back(symbol);
  store_initial_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_VAR].push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return ret;
}

std::string OptiNode::name_prefix() const {
  return "opti" + str(instance_number_) + "_";
}

Opti OptiNode::copy() const {
    return Opti::create(new OptiNode(*this));
}

void OptiNode::register_dual(MetaCon& c) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = "full";
  meta_data.n = c.canon.size1();
  meta_data.m = c.canon.size2();;
  meta_data.type = OPTI_DUAL_G;
  meta_data.count = count_++;
  meta_data.i = count_dual_++;

  MX symbol, ret;
  if (c.type==OPTI_PSD) {
    symbol = MX();
    ret = MX();
  } else {
    symbol = MX::sym(name_prefix()+"lam_g_"+str(count_dual_), c.canon.sparsity());

    casadi_assert_dev(c.canon.is_dense());

    Sparsity ret_sp = repmat(c.original.sparsity(), 1, c.n);

    casadi_int N = c.canon.sparsity().nnz();

    MX flat = vec(symbol);
    if (c.type==OPTI_DOUBLE_INEQUALITY) {
      MX v = MX::sym("v");
      MX decide_left_right = vertcat(if_else_zero(v<0, -v), if_else_zero(v>=0, v));
      Function sign = Function("sign", {v}, {decide_left_right});
      Function sign_map = sign.map(c.canon.sparsity().nnz());
      ret = MX(ret_sp, sign_map((c.flipped ? -1 : 1)*flat)[0].T());
    } else {
      casadi_int block_size = N / c.n;
      std::vector<MX> original_blocks = vertsplit(fabs(flat), block_size);
      std::vector<MX> blocks(N);
      for (casadi_int i=0;i<c.n;++i) {
        casadi_int p = c.flipped? c.n-i-1: i;
        blocks[p] = original_blocks[i];
      }
      ret = MX(ret_sp, vertcat(blocks));
    }
  }

  symbols_.push_back(symbol);
  store_initial_[OPTI_DUAL_G].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_DUAL_G].push_back(DM::nan(symbol.sparsity()));

  c.dual = ret;
  c.dual_canon = symbol;

  set_meta(symbol, meta_data);
}

MX OptiNode::parameter(casadi_int n, casadi_int m, const std::string& attribute) {
  casadi_assert_dev(attribute=="full");

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_PAR;
  meta_data.count = count_++;
  meta_data.i = count_par_++;

  MX symbol = MX::sym(name_prefix() + "p_" + str(count_par_), n, m);
  symbols_.push_back(symbol);
  store_initial_[OPTI_PAR].push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

Dict OptiNode::stats() const {
  assert_solved();
  return solver_.stats();
}

std::string OptiNode::return_status() const {
  Dict mystats;
  try {
    mystats = stats();
  } catch (...) {
    //
  }
  if (mystats.find("return_status")!=mystats.end())
    return mystats.at("return_status");
  return "unknown";
}

bool OptiNode::return_success() const {
  Dict mystats;
  try {
    mystats = stats();
  } catch (...) {
    //
  }
  if (mystats.find("success")!=mystats.end())
    return mystats.at("success");
  return false;
}

Function OptiNode::casadi_solver() const {
  return solver_;
}

void OptiNode::set_meta(const MX& m, const MetaVar& meta) {
  meta_[m.get()] = meta;
}

void OptiNode::set_meta_con(const MX& m, const MetaCon& meta) {
  meta_con_[m.get()] = meta;
}

void OptiNode::update_user_dict(const MX& m, const Dict& meta) {
  try {
    MetaCon m_update = get_meta_con(m);
    MetaVar m_update2 = get_meta(m_update.dual_canon);
    for (const auto & it : meta) {
        m_update.extra[it.first] = it.second;
        m_update2.extra[it.first] = it.second;
    }
    set_meta_con(m, m_update);
    set_meta(m_update.dual_canon, m_update2);
  } catch (exception& e) {
    for (const auto & s : MX::symvar(m)) {
      MetaVar m_update = get_meta(s);
      for (const auto & it : meta)
          m_update.extra[it.first] = it.second;
      set_meta(s, m_update);
    }
  }
}

Dict OptiNode::user_dict(const MX& m) const {
  try {
    MetaCon meta = get_meta_con(m);
    return meta.extra;
  } catch (exception& e) {
    MetaVar meta = get_meta(m);
    return meta.extra;
  }
}

MX OptiNode::dual(const MX& m) const {
  return meta_con(m).dual;
}

const MetaVar& OptiNode::meta(const MX& m) const {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

const MetaCon& OptiNode::meta_con(const MX& m) const {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

MetaVar& OptiNode::meta(const MX& m) {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

MetaCon& OptiNode::meta_con(const MX& m) {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

MetaVar OptiNode::get_meta(const MX& m) const {
  return meta(m);
}

MetaCon OptiNode::get_meta_con(const MX& m) const {
  return meta_con(m);
}

bool OptiNode::has(const MX& m) const {
  return meta_.find(m.get())!=meta_.end();
}

bool OptiNode::has_con(const MX& m) const {
  return meta_con_.find(m.get())!=meta_con_.end();
}

void OptiNode::assert_has(const MX& m) const {
  if (!has(m)) {
    VariableType vt;
    casadi_assert(m.is_symbolic(), "Symbol expected, got expression.");
    if (parse_opti_name(m.name(), vt)) {
      casadi_error("Unknown: " + describe(m));
    } else {
      casadi_error("Unknown: " + describe(m) + "\n"
        "Note: you cannot use a raw MX.sym in your Opti problem,"
        " only if you package it in a CasADi Function.");
    }
  }
}

void OptiNode::assert_has_con(const MX& m) const {
  casadi_assert(has_con(m), "Constraint not present in Opti stack.");
}

casadi_int OptiNode::instance_count_ = 0;

bool OptiNode::parse_opti_name(const std::string& name, VariableType& vt) const {
  casadi_int i = name.find("opti");
  if (i!=0) return false;

  i = name.find("_");
  i++;
  if (i==std::string::npos) return false;
  if (name.substr(i, 1)=="x") {
    vt = OPTI_VAR;
    return true;
  } else if (name.substr(i, 1)=="p") {
    vt = OPTI_PAR;
    return true;
  } else if (name.substr(i, 5)=="lam_g") {
    vt = OPTI_DUAL_G;
    return true;
  }

  return false;
}

std::string OptiNode::variable_type_to_string(VariableType vt) const {
  auto it = VariableType2String_.find(vt);
  if (it==VariableType2String_.end()) return "unknown variable type";
  return it->second;

}
std::map<VariableType, std::string> OptiNode::VariableType2String_ =
  {{OPTI_VAR, "decision variable"},
   {OPTI_PAR, "parameter"},
   {OPTI_DUAL_G, "dual variable"}};

std::vector<MX> OptiNode::initial() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR || meta(e).type==OPTI_DUAL_G)
      ret.push_back(e==store_initial_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiNode::value_variables() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR)
      ret.push_back(e==store_latest_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiNode::value_parameters() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_PAR)
      ret.push_back(e==store_initial_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

void OptiNode::bake() {
  casadi_assert(!f_.is_empty() || !g_.empty(),
    "You need to specify at least an objective (y calling 'minimize'), "
    "or a constraint (by calling 'subject_to').");

  symbol_active_.clear();
  symbol_active_.resize(symbols_.size());

  // Gather all expressions
  MX total_expr = vertcat(f_, veccat(g_));

  // Categorize the symbols appearing in those expressions
  for (const auto& d : symvar(total_expr))
    symbol_active_[meta(d).count] = true;

  std::vector<MX> x = active_symvar(OPTI_VAR);
  casadi_int offset = 0;
  for (const auto& v : x) {
    meta(v).start = offset;
    offset+= v.nnz();
    meta(v).stop = offset;
  }
  std::vector<MX> p = active_symvar(OPTI_PAR);

  // Fill the nlp definition
  nlp_["x"] = veccat(x);
  nlp_["p"] = veccat(p);

  nlp_["f"] = f_;

  offset = 0;
  for (casadi_int i=0;i<g_.size();++i) {
    MetaCon& r = meta_con(g_[i]);
    MetaVar& r2 = meta(r.dual_canon);
    symbol_active_[r2.count] = true;

    // Compute offsets for this constraint:
    // location into the global constraint variable
    r.start = offset;
    offset+= r.canon.nnz();
    r.stop = offset;

    r2.start = r.start;
    r2.stop  = r.stop;

  }

  lam_ = veccat(active_symvar(OPTI_DUAL_G));

  // Collect bounds and canonical form of constraints
  std::vector<MX> g_all;
  std::vector<MX> lbg_all;
  std::vector<MX> ubg_all;
  for (const auto& g : g_) {
    g_all.push_back(meta_con(g).canon);
    lbg_all.push_back(meta_con(g).lb);
    ubg_all.push_back(meta_con(g).ub);
  }

  nlp_["g"] = veccat(g_all);

  // Create bounds helper function
  MXDict bounds;
  bounds["p"] = nlp_["p"];
  bounds_lbg_ = veccat(lbg_all);
  bounds_ubg_ = veccat(ubg_all);

  bounds["lbg"] = bounds_lbg_;
  bounds["ubg"] = bounds_ubg_;

  bounds_ = Function("bounds", bounds, {"p"}, {"lbg", "ubg"});
  mark_problem_dirty(false);
}

void OptiNode::solver(const std::string& solver_name, const Dict& plugin_options,
                       const Dict& solver_options) {
  solver_name_ = solver_name;
  solver_options_ = plugin_options;
  if (!solver_options.empty())
    solver_options_[solver_name] = solver_options;
  mark_solver_dirty();
}

std::vector<MX> OptiNode::sort(const std::vector<MX>& v) const {
  // We exploit the fact that std::map is ordered

  // Populate map
  std::map<casadi_int, MX> unordered;
  for (const auto& d : v)
    unordered[meta(d).count] = d;

  // Loop over map (ordered)
  std::vector<MX> ret;
  for (auto const &e : unordered)
    ret.push_back(e.second);
  return ret;
}

std::vector<MX> OptiNode::symvar() const {
  return symbols_;
}

std::vector<MX> OptiNode::symvar(const MX& expr) const {
  return sort(MX::symvar(expr));
}

std::vector<MX> OptiNode::ineq_unchain(const MX& a, bool& flipped) {
  flipped = false;
  casadi_assert_dev(a.is_op(OP_LE) || a.is_op(OP_LT));

  // Is there inequalities in the left or right leaf?
  bool left  = a.dep(0).is_op(OP_LE) || a.dep(0).is_op(OP_LT);
  bool right = a.dep(1).is_op(OP_LE) || a.dep(1).is_op(OP_LT);
  casadi_assert_dev(!left || !right);

  if (!left && !right)
    return {a.dep(0), a.dep(1)}; // Simple inequality

  // We have a chain of inequalities
  bool ineq = !left;
  std::vector<MX> ret = {a.dep(!ineq)};
  MX e = a.dep(ineq);
  while (e.is_op(OP_LE) || e.is_op(OP_LT)) {
    casadi_assert_dev(!e.is_op(OP_EQ));
    casadi_assert_dev(!e.dep(!ineq).is_op(OP_LE) && !e.dep(!ineq).is_op(OP_LT));
    ret.push_back(e.dep(!ineq));
    e = e.dep(ineq);
  }
  ret.push_back(e);
  if (left) std::reverse(ret.begin(), ret.end());
  flipped = !left;

  return ret;
}

void OptiNode::assert_only_opti_symbols(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) assert_has(s);
}

void OptiNode::assert_only_opti_nondual(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) {
    assert_has(s);
    casadi_assert(meta(s).type!=OPTI_DUAL_G, "Dual variables forbidden in this context.");
  }
}

bool OptiNode::is_parametric(const MX& expr) const {
  return symvar(expr, OPTI_VAR).empty();
}

MetaCon OptiNode::canon_expr(const MX& expr) const {
  MX c = expr;

  MetaCon con;
  con.original = expr;

  if (c.is_op(OP_LE) || c.is_op(OP_LT)) { // Inequalities
    std::vector<MX> ret;
    bool flipped;
    std::vector<MX> args = ineq_unchain(c, flipped);
    std::vector<bool> parametric;
    for (auto &a : args) parametric.push_back(is_parametric(a));

    if (args.size()==2 && (parametric[0] || parametric[1])) {
      // case: g(x,p) <= bound(p)
      MX e = args[0]-args[1];
      if (e.is_vector()) {
        casadi_assert(!parametric[0] || !parametric[1],
          "Constraint must contain decision variables.");
        con.type = OPTI_INEQUALITY;
        if (parametric[0]) {
          con.lb = args[0]*DM::ones(e.sparsity());
          con.ub = inf*DM::ones(e.sparsity());
          con.canon = args[1]*DM::ones(e.sparsity());
        } else {
          con.lb = -inf*DM::ones(e.sparsity());
          con.ub = args[1]*DM::ones(e.sparsity());
          con.canon = args[0]*DM::ones(e.sparsity());
        }
        return con;
      }
      // Fall through to generic inequalities
    } else if (args.size()==3 && parametric[0] && parametric[2]) {
      // lb(p) <= g(x,p) <= ub(p)
      con.type = OPTI_DOUBLE_INEQUALITY;
      con.lb = args[0]*DM::ones(args[1].sparsity());
      con.ub = args[2]*DM::ones(args[1].sparsity());
      con.canon = args[1]*DM::ones(args[1].sparsity());
      con.flipped = flipped;
      con.n = 2;
      return con;
    }

    bool type_known = false;
    for (casadi_int j=0;j<args.size()-1;++j) {
      MX e = args[j]-args[j+1];
      if (e.is_vector()) {
        // g1(x,p) <= g2(x,p)
        ret.push_back(e);
        casadi_assert_dev(!type_known || con.type==OPTI_GENERIC_INEQUALITY);
        type_known = true;
        con.type = OPTI_GENERIC_INEQUALITY;
        con.flipped = flipped;
      } else {
        // A(x,p) >= b(p)
        e = args[j+1]-args[j];
        casadi_assert(e.size1()==e.size2(),
          "Matrix inequalities must be square. Did you mean element-wise inequality instead?");

        ret.push_back(e);
        casadi_assert_dev(!type_known || con.type==OPTI_PSD);
        type_known = true;
        con.type = OPTI_PSD;
      }
    }

    if (con.type==OPTI_GENERIC_INEQUALITY) {
      con.canon = veccat(ret);
      con.lb = -inf*DM::ones(con.canon.sparsity());
      con.ub = DM::zeros(con.canon.sparsity());
      con.n = ret.size();
    } else {
      con.canon = diagcat(ret);
      con.n = ret.size();
    }
    return con;
  } else if (c.is_op(OP_EQ)) { // Inequalities
    casadi_assert(!is_parametric(c.dep(0)) || !is_parametric(c.dep(1)),
      "Constraint must contain decision variables.");
    MX e = c.dep(0)-c.dep(1);
    if (is_parametric(c.dep(0))) {
      con.canon = c.dep(1)*DM::ones(e.sparsity());
      con.lb = c.dep(0)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
      casadi_assert(c.dep(0).size1()<=c.dep(1).size1() && c.dep(0).size2()<=c.dep(1).size2(),
        "Constraint shape mismatch.");
    } else if (is_parametric(c.dep(1))) {
      con.canon = c.dep(0)*DM::ones(e.sparsity());
      con.lb = c.dep(1)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
      casadi_assert(c.dep(1).size1()<=c.dep(0).size1() && c.dep(1).size2()<=c.dep(0).size2(),
        "Constraint shape mismatch.");
    } else {
      con.lb = DM::zeros(e.sparsity());
      con.canon = e;
      con.type = OPTI_GENERIC_EQUALITY;
    }
    con.ub = con.lb;
    return con;
  } else { // Something else
    con.type = OPTI_UNKNOWN;
    con.canon = c;
    return con;
  }

}

void OptiNode::assert_solved() const {
  casadi_assert(solved(),
    "This action is forbidden since you have not solved the Opti stack yet "
    "(with calling 'solve').");
}

void OptiNode::assert_baked() const {
  casadi_assert(!problem_dirty(),
    "This action is forbidden since you have not baked the Opti stack yet "
    "(with calling 'solve').");
}

void OptiNode::assert_empty() const {
  casadi_assert_dev(g_.empty());
  casadi_assert_dev(f_.is_empty());
}

void OptiNode::minimize(const MX& f) {
  assert_only_opti_nondual(f);
  mark_problem_dirty();
  casadi_assert(f.is_scalar(), "Objective must be scalar, got " + f.dim() + ".");
  f_ = f;
}

void OptiNode::subject_to(const MX& g) {
  assert_only_opti_nondual(g);
  mark_problem_dirty();
  g_.push_back(g);

  casadi_assert(!g.is_empty(),    "You passed an empty expression to `subject_to`. "
                                  "Make sure the number of rows and columns is non-zero. "
                                  "Got " + g.dim(true) + ".");
  casadi_assert(g.nnz()>0,        "You passed a fully sparse expression to `subject_to`. "
                                  "Make sure the expression has at least one nonzero. "
                                  "Got " + g.dim(true) + ".");
  casadi_assert(!g.is_constant(), "You passed a constant to `subject_to`. "
                                  "You need a symbol to form a constraint.");

  // Store the meta-data
  set_meta_con(g, canon_expr(g));
  register_dual(meta_con(g));
}

void OptiNode::subject_to() {
  mark_problem_dirty();
  g_.clear();
  store_initial_[OPTI_DUAL_G].clear();
  store_latest_[OPTI_DUAL_G].clear();
  count_dual_ = 0;
}

std::vector<MX> OptiNode::symvar(const MX& expr, VariableType type) const {
  std::vector<MX> ret;
  for (const auto& d : symvar(expr)) {
    if (meta(d).type==type) ret.push_back(d);
  }

  return ret;
}

void OptiNode::res(const DMDict& res) {
  const std::vector<double> & x_v = res.at("x").nonzeros();
  for (const auto &v : active_symvar(OPTI_VAR)) {
    casadi_int i = meta(v).i;
    std::vector<double> & data_v = store_latest_[OPTI_VAR][i].nonzeros();
    std::copy(x_v.begin()+meta(v).start, x_v.begin()+meta(v).stop, data_v.begin());
  }
  if (res.find("lam_g")!=res.end()) {
    const std::vector<double> & lam_v = res.at("lam_g").nonzeros();
    for (const auto &v : active_symvar(OPTI_DUAL_G)) {
      casadi_int i = meta(v).i;
      std::vector<double> & data_v = store_latest_[OPTI_DUAL_G][i].nonzeros();
      std::copy(lam_v.begin()+meta(v).start, lam_v.begin()+meta(v).stop, data_v.begin());
    }
  }
  res_ = res;
  mark_solved();
}

bool OptiNode::old_callback() const {
  if (callback_.is_null()) return false;
  InternalOptiCallback* cb = static_cast<InternalOptiCallback*>(callback_.get());
  return !cb->associated_with(this);
}
// Solve the problem
OptiSol OptiNode::solve() {

  if (problem_dirty()) {
    bake();
  }
  // Verify the constraint types
  for (const auto& g : g_) {
    if (meta_con(g).type==OPTI_PSD)
      casadi_error("Psd constraints not implemented yet. "
      "Perhaps you intended an element-wise inequality? "
      "In that case, make sure that the matrix is flattened (e.g. mat(:)).");
  }

  bool solver_update =  solver_dirty() || old_callback() || (user_callback_ && callback_.is_null());

  if (solver_update) {
    Dict opts = solver_options_;

    // Handle callbacks
    if (user_callback_) {
      callback_ = Function::create(new InternalOptiCallback(*this), Dict());
      opts["iteration_callback"] = callback_;
    }

    casadi_assert(solver_name_!="",
      "You must call 'solver' on the Opti stack to select a solver. "
      "Suggestion: opti.solver('ipopt')");

    solver_ = nlpsol("solver", solver_name_, nlp_, opts);
    mark_solver_dirty(false);
  }

  solve_prepare();
  res(solve_actual(arg_));

  std::string ret = return_status();

  casadi_assert(return_success(),
    "Solver failed. You may use opti.debug.value to investigate the latest values of variables."
    " return_status is '" + ret + "'");

  return Opti(this);
}

// Solve the problem
void OptiNode::solve_prepare() {


  // Verify the constraint types
  for (const auto& g : g_) {
    if (meta_con(g).type==OPTI_UNKNOWN)
     casadi_error("Constraint type unknown. Use ==, >= or <= .");
  }

  if (user_callback_) {
    InternalOptiCallback* cb = static_cast<InternalOptiCallback*>(callback_.get());
    cb->reset();
  }

  // Get initial guess and parameter values
  arg_["x0"]     = veccat(active_values(OPTI_VAR));
  arg_["p"]      = veccat(active_values(OPTI_PAR));
  arg_["lam_g0"] = veccat(active_values(OPTI_DUAL_G));
  if (!arg_["p"].is_regular()) {
    std::vector<MX> s = active_symvar(OPTI_PAR);
    std::vector<DM> v = active_values(OPTI_PAR);
    for (casadi_int i=0;i<s.size();++i) {
      casadi_assert(v[i].is_regular(),
        "You have forgotten to assign a value to a parameter ('set_value'), "
        "or have set it to NaN/Inf:\n" + describe(s[i], 1));
    }
  }

  // Evaluate bounds for given parameter values
  DMDict arg;
  arg["p"] = arg_["p"];
  DMDict res = bounds_(arg);
  arg_["lbg"] = res["lbg"];
  arg_["ubg"] = res["ubg"];

}

DMDict OptiNode::solve_actual(const DMDict& arg) {
  return solver_(arg);
}

bool override_num(const std::map<casadi_int, MX> & temp, std::vector<DM>& num, casadi_int i) {
  // Override when values are supplied
  auto it = temp.find(i);
  if (it==temp.end()) {
    return true;
  } else {
    Slice all;
    DM t = static_cast<DM>(it->second);
    num.back().set(t, false, all, all);
  }
  return false;
}

DM OptiNode::value(const MX& expr, const std::vector<MX>& values) const {
  std::vector<MX> x   = symvar(expr, OPTI_VAR);
  std::vector<MX> p   = symvar(expr, OPTI_PAR);
  std::vector<MX> lam = symvar(expr, OPTI_DUAL_G);

  Function helper = Function("helper", std::vector<MX>{veccat(x), veccat(p), veccat(lam)}, {expr});
  if (helper.has_free())
    casadi_error("This expression has symbols that are not defined "
      "within Opti using variable/parameter.");

  std::map<VariableType, std::map<casadi_int, MX> > temp;
  temp[OPTI_DUAL_G] = std::map<casadi_int, MX>();
  for (const auto& v : values) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_int i = meta(v.dep(1)).i;
    casadi_assert_dev(v.dep(0).is_constant());
    temp[meta(v.dep(1)).type][i] = v.dep(0);
  }

  bool undecided_vars = false;
  std::vector<DM> x_num;
  for (const auto& e : x) {
    casadi_int i = meta(e).i;
    x_num.push_back(store_latest_.at(OPTI_VAR).at(i));
    undecided_vars |= override_num(temp[OPTI_VAR], x_num, i);
  }

  std::vector<DM> lam_num;
  for (const auto& e : lam) {
    casadi_int i = meta(e).i;
    casadi_assert(i<store_latest_.at(OPTI_DUAL_G).size(),
      "This expression has a dual for a constraint that is not given to Opti:\n" +
      describe(e, 1));
    lam_num.push_back(store_latest_.at(OPTI_DUAL_G).at(i));
    undecided_vars |= override_num(temp[OPTI_DUAL_G], lam_num, i);
  }

  std::vector<DM> p_num;
  for (const auto& e : p) {
    casadi_int i = meta(e).i;
    p_num.push_back(store_initial_.at(OPTI_PAR).at(i));
    override_num(temp[OPTI_PAR], p_num, i);
    casadi_assert(p_num.back().is_regular(),
      "This expression depends on a parameter with unset value:\n"+
      describe(e, 1));
  }

  if (undecided_vars) {
    assert_solved();
    for (const auto& e : x)
      casadi_assert(symbol_active_[meta(e).count],
        "This expression has symbols that do not appear in the constraints and objective:\n" +
        describe(e, 1));
    for (const auto& e : lam)
      casadi_assert(symbol_active_[meta(e).count],
        "This expression has a dual for a constraint that is not given to Opti:\n" +
        describe(e, 1));
  }

  std::vector<DM> arg = helper(std::vector<DM>{veccat(x_num), veccat(p_num), veccat(lam_num)});
  return arg[0];
}

void OptiNode::assert_active_symbol(const MX& m) const {
  assert_has(m);
  assert_baked();
  casadi_assert(symbol_active_[meta(m).count], "Opti symbol is not used in Solver."
    " It does not make sense to assign a value to it:\n" + describe(m, 1));
}

void OptiNode::set_initial(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_assert_dev(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_initial(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiNode::set_value(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_assert_dev(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_value(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiNode::set_value_internal(const MX& x, const DM& v) {
  mark_solved(false);
  casadi_assert_dev(v.is_regular());
  if (x.is_symbolic()) {
    DM& target = store_initial_[meta(x).type][meta(x).i];
    Slice all;
    target.set(v, false, all, all);
    return;
  }

  // Obtain symbolic primitives
  std::vector<MX> symbols = MX::symvar(x);
  MX symbols_cat = veccat(symbols);

  std::string failmessage = "You cannot set initial/value of an arbitrary expression. "
    "Use symbols or simple mappings of symbols.";

  // Assert x is linear in its symbolic primitives
  for (bool b : which_depends(x, symbols_cat, 2, false)) casadi_assert(!b, failmessage);

  // Evaluate jacobian of expr wrt symbols
  Dict opts = {{"compact", true}};
  Function Jf("Jf", std::vector<MX>{}, std::vector<MX>{jacobian(x, veccat(symbols), opts)});
  DM J = Jf(std::vector<DM>{})[0];
  Sparsity sp_JT = J.T().sparsity();

  Function Ff("Ff", symbols, {x});
  DM E = Ff(std::vector<DM>(symbols.size(), 0))[0];
  std::vector<double>& e = E.nonzeros();

  // Cast the v input into the expected sparsity
  Slice all;
  DM value(x.sparsity());
  value.set(v, false, all, all);

  // Purge empty rows
  std::vector<casadi_int> filled_rows = sum2(J).get_row();
  J = J(filled_rows, all);

  // Get rows and columns of the mapping
  std::vector<casadi_int> row, col;
  J.sparsity().get_triplet(row, col);
  const std::vector<double>& scaling = J.nonzeros();
  const std::vector<double>& data_original = value.nonzeros();

  std::vector<double> data; data.reserve(value.nnz());
  for (casadi_int i=0;i<value.nnz();++i) {
    double v = data_original[i];
    casadi_int nz = sp_JT.colind()[i+1]-sp_JT.colind()[i];
    casadi_assert(nz<=1, failmessage);
    if (nz) {
      data.push_back(v);
    } else {
      casadi_assert(v==e[i], "In initial/value assignment: "
        "inconsistent numerical values. At nonzero " + str(i) + ", lhs has "
        + str(e[i]) + ", while rhs has " + str(v) + ".");
    }
  }

  // Contiguous workspace for nonzeros of all involved symbols
  std::vector<double> temp(symbols_cat.nnz(), casadi::nan);
  for (casadi_int k=0;k<data.size();++k) {
    double& lhs = temp[col[k]];
    double rhs = data[row[k]]/scaling[row[k]];
    if (std::isnan(lhs)) {
      // Assign in the workspace
      lhs = rhs;
    } else {
      casadi_assert(lhs==rhs, "Initial/value assignment with mapping is ambiguous.");
    }
  }

  casadi_int offset = 0;
  for (const auto & s : symbols) {
    DM& target = store_initial_[meta(s).type][meta(s).i];
    std::vector<double>& data = target.nonzeros();
    // Loop over nonzeros in each symbol
    for (casadi_int i=0;i<s.nnz();++i) {
      // Copy from the workspace (barring fields that were not set)
      double v = temp[offset+i];
      if (!std::isnan(v)) data[i] = v;
    }
    offset+=s.nnz();
  }

}

void OptiNode::set_initial(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x))
    casadi_assert(meta(s).type!=OPTI_PAR,
      "You cannot set an initial value for a parameter. Did you mean 'set_value'?");
  set_value_internal(x, v);
}

void OptiNode::set_value(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x))
    casadi_assert(meta(s).type!=OPTI_VAR,
      "You cannot set a value for a variable. Did you mean 'set_initial'?");
  set_value_internal(x, v);
}

std::vector<MX> OptiNode::active_symvar(VariableType type) const {
  if (symbol_active_.empty()) return std::vector<MX>{};
  std::vector<MX> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type)
      ret.push_back(s);
  }
  return ret;
}

std::vector<DM> OptiNode::active_values(VariableType type) const {
  if (symbol_active_.empty()) return std::vector<DM>{};
  std::vector<DM> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type) {
      ret.push_back(store_initial_.at(meta(s).type)[meta(s).i]);
    }
  }
  return ret;
}

void OptiNode::disp(ostream &stream, bool more) const {

}

casadi_int OptiNode::instance_number() const {
    return instance_number_;
}

} // namespace casadi
