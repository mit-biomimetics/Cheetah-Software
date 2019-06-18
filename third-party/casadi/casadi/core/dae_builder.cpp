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


#include "dae_builder.hpp"

#include <map>
#include <set>
#include <string>
#include <sstream>
#include <ctime>
#include <cctype>

#include "casadi_misc.hpp"
#include "exception.hpp"
#include "code_generator.hpp"
#include "calculus.hpp"
#include "xml_file.hpp"
#include "external.hpp"

using namespace std;
namespace casadi {

  DaeBuilder::DaeBuilder() {
    this->t = MX::sym("t");
  }

  void DaeBuilder::parse_fmi(const std::string& filename) {

    // Load
    XmlFile xml_file("tinyxml");
    XmlNode document = xml_file.parse(filename);

    // **** Add model variables ****
    {
      // Get a reference to the ModelVariables node
      const XmlNode& modvars = document[0]["ModelVariables"];

      // Add variables
      for (casadi_int i=0; i<modvars.size(); ++i) {

        // Get a reference to the variable
        const XmlNode& vnode = modvars[i];

        // Get the attributes
        string name        = vnode.getAttribute("name");
        casadi_int valueReference;
        vnode.readAttribute("valueReference", valueReference);
        string variability = vnode.getAttribute("variability");
        string causality   = vnode.getAttribute("causality");
        string alias       = vnode.getAttribute("alias");

        // Skip to the next variable if its an alias
        if (alias.compare("alias") == 0 || alias.compare("negatedAlias") == 0)
          continue;

        // Get the name
        const XmlNode& nn = vnode["QualifiedName"];
        string qn = qualified_name(nn);

        // Add variable, if not already added
        if (varmap_.find(qn)==varmap_.end()) {

          // Create variable
          Variable var(name);

          // Value reference
          var.valueReference = valueReference;

          // Variability
          if (variability.compare("constant")==0)
            var.variability = CONSTANT;
          else if (variability.compare("parameter")==0)
            var.variability = PARAMETER;
          else if (variability.compare("discrete")==0)
            var.variability = DISCRETE;
          else if (variability.compare("continuous")==0)
            var.variability = CONTINUOUS;
          else
            throw CasadiException("Unknown variability");

          // Causality
          if (causality.compare("input")==0)
            var.causality = INPUT;
          else if (causality.compare("output")==0)
            var.causality = OUTPUT;
          else if (causality.compare("internal")==0)
            var.causality = INTERNAL;
          else
            throw CasadiException("Unknown causality");

          // Alias
          if (alias.compare("noAlias")==0)
            var.alias = NO_ALIAS;
          else if (alias.compare("alias")==0)
            var.alias = ALIAS;
          else if (alias.compare("negatedAlias")==0)
            var.alias = NEGATED_ALIAS;
          else
            throw CasadiException("Unknown alias");

          // Other properties
          if (vnode.hasChild("Real")) {
            const XmlNode& props = vnode["Real"];
            props.readAttribute("unit", var.unit, false);
            props.readAttribute("displayUnit", var.display_unit, false);
            props.readAttribute("min", var.min, false);
            props.readAttribute("max", var.max, false);
            props.readAttribute("initialGuess", var.guess, false);
            props.readAttribute("start", var.start, false);
            props.readAttribute("nominal", var.nominal, false);
            props.readAttribute("free", var.free, false);
          }

          // Variable category
          if (vnode.hasChild("VariableCategory")) {
            string cat = vnode["VariableCategory"].getText();
            if (cat.compare("derivative")==0)
              var.category = CAT_DERIVATIVE;
            else if (cat.compare("state")==0)
              var.category = CAT_STATE;
            else if (cat.compare("dependentConstant")==0)
              var.category = CAT_DEPENDENT_CONSTANT;
            else if (cat.compare("independentConstant")==0)
              var.category = CAT_INDEPENDENT_CONSTANT;
            else if (cat.compare("dependentParameter")==0)
              var.category = CAT_DEPENDENT_PARAMETER;
            else if (cat.compare("independentParameter")==0)
              var.category = CAT_INDEPENDENT_PARAMETER;
            else if (cat.compare("algebraic")==0)
              var.category = CAT_ALGEBRAIC;
            else
              throw CasadiException("Unknown variable category: " + cat);
          }

          // Add to list of variables
          add_variable(qn, var);

          // Sort expression
          switch (var.category) {
          case CAT_DERIVATIVE:
            // Skip - meta information about time derivatives is
            //        kept together with its parent variable
            break;
          case CAT_STATE:
            this->s.push_back(var.v);
            this->sdot.push_back(var.d);
            break;
          case CAT_DEPENDENT_CONSTANT:
            // Skip
            break;
          case CAT_INDEPENDENT_CONSTANT:
            // Skip
            break;
          case CAT_DEPENDENT_PARAMETER:
            // Skip
            break;
          case CAT_INDEPENDENT_PARAMETER:
            if (var.free) {
              this->p.push_back(var.v);
            } else {
              // Skip
            }
            break;
          case CAT_ALGEBRAIC:
            if (var.causality == INTERNAL) {
              this->s.push_back(var.v);
              this->sdot.push_back(var.d);
            } else if (var.causality == INPUT) {
              this->u.push_back(var.v);
            }
            break;
          default:
            casadi_error("Unknown category");
          }
        }
      }
    }

    // **** Add binding equations ****
    {
      // Get a reference to the BindingEquations node
      const XmlNode& bindeqs = document[0]["equ:BindingEquations"];

      for (casadi_int i=0; i<bindeqs.size(); ++i) {
        const XmlNode& beq = bindeqs[i];

        // Get the variable and binding expression
        Variable& var = read_variable(beq[0]);
        MX bexpr = read_expr(beq[1][0]);
        this->d.push_back(var.v);
        this->ddef.push_back(bexpr);
      }
    }

    // **** Add dynamic equations ****
    {
      // Get a reference to the DynamicEquations node
      const XmlNode& dyneqs = document[0]["equ:DynamicEquations"];

      // Add equations
      for (casadi_int i=0; i<dyneqs.size(); ++i) {

        // Get a reference to the variable
        const XmlNode& dnode = dyneqs[i];

        // Add the differential equation
        MX de_new = read_expr(dnode[0]);
        this->dae.push_back(de_new);
      }
    }

    // **** Add initial equations ****
    {
      // Get a reference to the DynamicEquations node
      const XmlNode& initeqs = document[0]["equ:InitialEquations"];

      // Add equations
      for (casadi_int i=0; i<initeqs.size(); ++i) {

        // Get a reference to the node
        const XmlNode& inode = initeqs[i];

        // Add the differential equations
        for (casadi_int i=0; i<inode.size(); ++i) {
          this->init.push_back(read_expr(inode[i]));
        }
      }
    }

    // **** Add optimization ****
    if (document[0].hasChild("opt:Optimization")) {

      // Get a reference to the DynamicEquations node
      const XmlNode& opts = document[0]["opt:Optimization"];
      for (casadi_int i=0; i<opts.size(); ++i) {

        // Get a reference to the node
        const XmlNode& onode = opts[i];

        // Get the type
        if (onode.checkName("opt:ObjectiveFunction")) { // mayer term
          try {
            // Add components
            for (casadi_int i=0; i<onode.size(); ++i) {
              const XmlNode& var = onode[i];

              // If string literal, ignore
              if (var.checkName("exp:StringLiteral"))
                continue;

              // Read expression
              MX v = read_expr(var);

              // Treat as an output
              add_y("mterm", v);
            }
          } catch(exception& ex) {
            throw CasadiException(std::string("addObjectiveFunction failed: ") + ex.what());
          }
        } else if (onode.checkName("opt:IntegrandObjectiveFunction")) {
          try {
            for (casadi_int i=0; i<onode.size(); ++i) {
              const XmlNode& var = onode[i];

              // If string literal, ignore
              if (var.checkName("exp:StringLiteral")) continue;

              // Read expression
              MX v = read_expr(var);

              // Treat as a quadrature state
              add_q("lterm");
              add_quad("lterm_rhs", v);
            }
          } catch(exception& ex) {
            throw CasadiException(std::string("addIntegrandObjectiveFunction failed: ")
                                  + ex.what());
          }
        } else if (onode.checkName("opt:IntervalStartTime")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:IntervalFinalTime")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:TimePoints")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:PointConstraints")) {
          casadi_warning("opt:PointConstraints not supported, ignored");
        } else if (onode.checkName("opt:Constraints")) {
          casadi_warning("opt:Constraints not supported, ignored");
        } else if (onode.checkName("opt:PathConstraints")) {
          casadi_warning("opt:PointConstraints not supported, ignored");
        } else {
          casadi_warning("DaeBuilder::addOptimization: Unknown node " + str(onode.name()));
        }
      }
    }

    // Make sure that the dimensions are consistent at this point
    if (this->s.size()!=this->dae.size()) {
      casadi_warning("The number of differential-algebraic equations does not match "
                     "the number of implicitly defined states.");
    }
    if (this->z.size()!=this->alg.size()) {
      casadi_warning("The number of algebraic equations (equations not involving "
                    "differentiated variables) does not match the number of "
                    "algebraic variables.");
    }
  }

  Variable& DaeBuilder::read_variable(const XmlNode& node) {
    // Qualified name
    string qn = qualified_name(node);

    // Find and return the variable
    return variable(qn);
  }

  MX DaeBuilder::read_expr(const XmlNode& node) {
    const string& fullname = node.name();
    if (fullname.find("exp:")== string::npos) {
      casadi_error("DaeBuilder::read_expr: unknown - expression is supposed to "
                   "start with 'exp:' , got " + fullname);
    }

    // Chop the 'exp:'
    string name = fullname.substr(4);

    // The switch below is alphabetical, and can be thus made more efficient,
    // for example by using a switch statement of the first three letters,
    // if it would ever become a bottleneck
    if (name.compare("Add")==0) {
      return read_expr(node[0]) + read_expr(node[1]);
    } else if (name.compare("Acos")==0) {
      return acos(read_expr(node[0]));
    } else if (name.compare("Asin")==0) {
      return asin(read_expr(node[0]));
    } else if (name.compare("Atan")==0) {
      return atan(read_expr(node[0]));
    } else if (name.compare("Cos")==0) {
      return cos(read_expr(node[0]));
    } else if (name.compare("Der")==0) {
      const Variable& v = read_variable(node[0]);
      return v.d;
    } else if (name.compare("Div")==0) {
      return read_expr(node[0]) / read_expr(node[1]);
    } else if (name.compare("Exp")==0) {
      return exp(read_expr(node[0]));
    } else if (name.compare("Identifier")==0) {
      return read_variable(node).v;
    } else if (name.compare("IntegerLiteral")==0) {
      casadi_int val;
      node.getText(val);
      return val;
    } else if (name.compare("Instant")==0) {
      double val;
      node.getText(val);
      return val;
    } else if (name.compare("Log")==0) {
      return log(read_expr(node[0]));
    } else if (name.compare("LogLeq")==0) { // Logical less than equal
      return read_expr(node[0]) <= read_expr(node[1]);
    } else if (name.compare("LogGeq")==0) { // Logical greater than equal
      return read_expr(node[0]) >= read_expr(node[1]);
    } else if (name.compare("LogLt")==0) { // Logical less than
      return read_expr(node[0]) < read_expr(node[1]);
    } else if (name.compare("LogGt")==0) { // Logical greater than
      return read_expr(node[0]) > read_expr(node[1]);
    } else if (name.compare("Max")==0) {
      return fmax(read_expr(node[0]), read_expr(node[1]));
    } else if (name.compare("Min")==0) {
      return fmin(read_expr(node[0]), read_expr(node[1]));
    } else if (name.compare("Mul")==0) { // Multiplication
      return read_expr(node[0]) * read_expr(node[1]);
    } else if (name.compare("Neg")==0) {
      return -read_expr(node[0]);
    } else if (name.compare("NoEvent")==0) {
      // NOTE: This is a workaround, we assume that whenever NoEvent occurs,
      // what is meant is a switch
      casadi_int n = node.size();

      // Default-expression
      MX ex = read_expr(node[n-1]);

      // Evaluate ifs
      for (casadi_int i=n-3; i>=0; i -= 2) {
        ex = if_else(read_expr(node[i]), read_expr(node[i+1]), ex);
      }

      return ex;
    } else if (name.compare("Pow")==0) {
      return pow(read_expr(node[0]), read_expr(node[1]));
    } else if (name.compare("RealLiteral")==0) {
      double val;
      node.getText(val);
      return val;
    } else if (name.compare("Sin")==0) {
      return sin(read_expr(node[0]));
    } else if (name.compare("Sqrt")==0) {
      return sqrt(read_expr(node[0]));
    } else if (name.compare("StringLiteral")==0) {
      throw CasadiException(node.getText());
    } else if (name.compare("Sub")==0) {
      return read_expr(node[0]) - read_expr(node[1]);
    } else if (name.compare("Tan")==0) {
      return tan(read_expr(node[0]));
    } else if (name.compare("Time")==0) {
      return t;
    } else if (name.compare("TimedVariable")==0) {
      return read_variable(node[0]).v;
    }

    // throw error if reached this point
    throw CasadiException(string("DaeBuilder::read_expr: Unknown node: ") + name);

  }

  void DaeBuilder::disp(std::ostream& stream, bool more) const {
    // Assert correctness
    if (more) sanity_check();

    // Print dimensions
    stream << "ns = " << this->s.size() << ", "
           << "nx = " << this->x.size() << ", "
           << "nz = " << this->z.size() << ", "
           << "nq = " << this->q.size() << ", "
           << "ny = " << this->y.size() << ", "
           << "np = " << this->p.size() << ", "
           << "nd = " << this->d.size() << ", "
           << "nu = " << this->u.size();

    // Quick return?
    if (!more) return;
    stream << endl;

    // Print the functions
    if (!fun_.empty()) {
      stream << "Functions" << endl;
      for (const Function& f : fun_) {
        stream << "  " << f << endl;
      }
    }

    // Print the variables
    stream << "Variables" << endl;
    stream << "  t = " << str(this->t) << endl;
    if (!this->s.empty()) stream << "  s = " << str(this->s) << endl;
    if (!this->x.empty()) stream << "  x = " << str(this->x) << endl;
    if (!this->z.empty()) stream << "  z =  " << str(this->z) << endl;
    if (!this->q.empty()) stream << "  q =  " << str(this->q) << endl;
    if (!this->y.empty()) stream << "  y =  " << str(this->y) << endl;
    if (!this->p.empty()) stream << "  p =  " << str(this->p) << endl;
    if (!this->d.empty()) stream << "  d =  " << str(this->d) << endl;
    if (!this->u.empty()) stream << "  u =  " << str(this->u) << endl;

    if (!this->d.empty()) {
      stream << "Dependent parameters" << endl;
      for (casadi_int i=0; i<this->d.size(); ++i)
        stream << "  " << str(this->d[i]) << " == " << str(this->ddef[i]) << endl;
    }

    if (!this->dae.empty()) {
      stream << "Fully-implicit differential-algebraic equations" << endl;
      for (casadi_int k=0; k<this->dae.size(); ++k) {
        stream << "  0 == " << this->dae[k] << endl;
      }
    }

    if (!this->x.empty()) {
      stream << "Differential equations" << endl;
      for (casadi_int k=0; k<this->x.size(); ++k) {
        stream << "  " << str(der(this->x[k])) << " == " << str(this->ode[k]) << endl;
      }
    }

    if (!this->alg.empty()) {
      stream << "Algebraic equations" << endl;
      for (casadi_int k=0; k<this->z.size(); ++k) {
        stream << "  0 == " << str(this->alg[k]) << endl;
      }
    }

    if (!this->q.empty()) {
      stream << "Quadrature equations" << endl;
      for (casadi_int k=0; k<this->q.size(); ++k) {
        stream << "  " << str(der(this->q[k])) << " == " << str(this->quad[k]) << endl;
      }
    }

    if (!this->init.empty()) {
      stream << "Initial equations" << endl;
      for (casadi_int k=0; k<this->init.size(); ++k) {
        stream << "  0 == " << str(this->init[k]) << endl;
      }
    }

    if (!this->y.empty()) {
      stream << "Output variables" << endl;
      for (casadi_int i=0; i<this->y.size(); ++i) {
        stream << "  " << str(this->y[i]) << " == " << str(this->ydef[i]) << endl;
      }
    }
  }

  void DaeBuilder::eliminate_quad() {
    // Move all the quadratures to the list of differential states
    this->x.insert(this->x.end(), this->q.begin(), this->q.end());
    this->q.clear();
  }

  void DaeBuilder::scale_variables() {
    // Assert correctness
    sanity_check();

    // Gather variables and expressions to replace
    vector<MX> v_id, v_rep;
    for (VarMap::iterator it=varmap_.begin(); it!=varmap_.end(); ++it) {
      if (it->second.nominal!=1) {
        Variable& v=it->second;
        casadi_assert_dev(v.nominal!=0);
        v.min /= v.nominal;
        v.max /= v.nominal;
        v.start /= v.nominal;
        v.derivative_start /= v.nominal;
        v.guess /= v.nominal;
        v_id.push_back(v.v);
        v_id.push_back(v.d);
        v_rep.push_back(v.v * v.nominal);
        v_rep.push_back(v.d * v.nominal);
      }
    }

    // Quick return if no expressions to substitute
    if (v_id.empty()) return;

    // Collect all expressions to be replaced
    vector<MX> ex;
    ex.insert(ex.end(), this->ode.begin(), this->ode.end());
    ex.insert(ex.end(), this->dae.begin(), this->dae.end());
    ex.insert(ex.end(), this->alg.begin(), this->alg.end());
    ex.insert(ex.end(), this->quad.begin(), this->quad.end());
    ex.insert(ex.end(), this->ddef.begin(), this->ddef.end());
    ex.insert(ex.end(), this->ydef.begin(), this->ydef.end());
    ex.insert(ex.end(), this->init.begin(), this->init.end());

    // Substitute all at once (more efficient since they may have common subexpressions)
    ex = substitute(ex, v_id, v_rep);

    // Get the modified expressions
    vector<MX>::const_iterator it=ex.begin();
    for (casadi_int i=0; i<this->x.size(); ++i) this->ode[i] = *it++ / nominal(this->x[i]);
    for (casadi_int i=0; i<this->s.size(); ++i) this->dae[i] = *it++;
    for (casadi_int i=0; i<this->z.size(); ++i) this->alg[i] = *it++;
    for (casadi_int i=0; i<this->q.size(); ++i) this->quad[i] = *it++ / nominal(this->q[i]);
    for (casadi_int i=0; i<this->d.size(); ++i) this->ddef[i] = *it++ / nominal(this->d[i]);
    for (casadi_int i=0; i<this->y.size(); ++i) this->ydef[i] = *it++ / nominal(this->y[i]);
    for (casadi_int i=0; i<this->init.size(); ++i) this->init[i] = *it++;
    casadi_assert_dev(it==ex.end());

    // Nominal value is 1 after scaling
    for (VarMap::iterator it=varmap_.begin(); it!=varmap_.end(); ++it) {
      it->second.nominal=1;
    }
  }

  void DaeBuilder::sort_d() {
    // Quick return if no intermediates
    if (this->d.empty()) return;

    // Find out which intermediates depends on which other
    Function f("tmp", {vertcat(this->d)}, {vertcat(this->d) - vertcat(this->ddef)});
    Sparsity sp = f.sparsity_jac(0, 0);
    casadi_assert_dev(sp.is_square());

    // BLT transformation
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Resort equations and variables
    vector<MX> ddefnew(this->d.size()), dnew(this->d.size());
    for (casadi_int i=0; i<colperm.size(); ++i) {
      // Permute equations
      ddefnew[i] = this->ddef[colperm[i]];

      // Permute variables
      dnew[i] = this->d[colperm[i]];
    }
    this->ddef = ddefnew;
    this->d = dnew;
  }

  void DaeBuilder::split_d() {
    // Quick return if no intermediates
    if (this->d.empty()) return;

    // Begin by sorting the dependent parameters
    sort_d();

    // Sort the equations by causality
    vector<MX> ex;
    substitute_inplace(this->d, this->ddef, ex);

    // Make sure that the interdependencies have been properly eliminated
    casadi_assert_dev(!depends_on(vertcat(this->ddef), vertcat(this->d)));
  }

  void DaeBuilder::eliminate_d() {
    // Quick return if possible
    if (this->d.empty()) return;

    // Begin by sorting the dependent parameters
    sort_d();

    // Collect all expressions to be replaced
    vector<MX> ex;
    ex.insert(ex.end(), this->ode.begin(), this->ode.end());
    ex.insert(ex.end(), this->dae.begin(), this->dae.end());
    ex.insert(ex.end(), this->alg.begin(), this->alg.end());
    ex.insert(ex.end(), this->quad.begin(), this->quad.end());
    ex.insert(ex.end(), this->ydef.begin(), this->ydef.end());
    ex.insert(ex.end(), this->init.begin(), this->init.end());

    // Substitute all at once (since they may have common subexpressions)
    substitute_inplace(this->d, this->ddef, ex);

    // Get the modified expressions
    vector<MX>::const_iterator it=ex.begin();
    for (casadi_int i=0; i<this->x.size(); ++i) this->ode[i] = *it++;
    for (casadi_int i=0; i<this->s.size(); ++i) this->dae[i] = *it++;
    for (casadi_int i=0; i<this->z.size(); ++i) this->alg[i] = *it++;
    for (casadi_int i=0; i<this->q.size(); ++i) this->quad[i] = *it++;
    for (casadi_int i=0; i<this->y.size(); ++i) this->ydef[i] = *it++;
    for (casadi_int i=0; i<this->init.size(); ++i) this->init[i] = *it++;
    casadi_assert_dev(it==ex.end());
  }

  void DaeBuilder::scale_equations() {
    casadi_error("DaeBuilder::scale_equations broken");
  }

  void DaeBuilder::sort_dae() {
    // Quick return if no differential states
    if (this->x.empty()) return;

    // Find out which differential equation depends on which differential state
    Function f("tmp", {vertcat(this->sdot)}, {vertcat(this->dae)});
    Sparsity sp = f.sparsity_jac(0, 0);
    casadi_assert_dev(sp.is_square());

    // BLT transformation
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Resort equations and variables
    vector<MX> daenew(this->s.size()), snew(this->s.size()), sdotnew(this->s.size());
    for (casadi_int i=0; i<rowperm.size(); ++i) {
      // Permute equations
      daenew[i] = this->dae[rowperm[i]];

      // Permute variables
      snew[i] = this->s[colperm[i]];
      sdotnew[i] = this->sdot[colperm[i]];
    }
    this->dae = daenew;
    this->s = snew;
    this->sdot = sdotnew;
  }

  void DaeBuilder::sort_alg() {
    // Quick return if no algebraic states
    if (this->z.empty()) return;

    // Find out which algebraic equation depends on which algebraic state
    Function f("tmp", {vertcat(this->z)}, {vertcat(this->alg)});
    Sparsity sp = f.sparsity_jac(0, 0);
    casadi_assert_dev(sp.is_square());

    // BLT transformation
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Resort equations and variables
    vector<MX> algnew(this->z.size()), znew(this->z.size());
    for (casadi_int i=0; i<rowperm.size(); ++i) {
      // Permute equations
      algnew[i] = this->alg[rowperm[i]];

      // Permute variables
      znew[i] = this->z[colperm[i]];
    }
    this->alg = algnew;
    this->z = znew;
  }

  void DaeBuilder::make_semi_explicit() {
    // Only works if there are no i
    eliminate_d();

    // Separate the algebraic variables and equations
    split_dae();

    // Quick return if there are no implicitly defined states
    if (this->s.empty()) return;

    // Write the ODE as a function of the state derivatives
    Function f("tmp", {vertcat(this->sdot)}, {vertcat(this->dae)});

    // Get the sparsity of the Jacobian which can be used to determine which
    // variable can be calculated from which other
    Sparsity sp = f.sparsity_jac(0, 0);
    casadi_assert_dev(sp.is_square());

    // BLT transformation
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    casadi_int nb = sp.btf(rowperm, colperm, rowblock, colblock,
                                  coarse_rowblock, coarse_colblock);

    // Resort equations and variables
    vector<MX> daenew(this->s.size()), snew(this->s.size()), sdotnew(this->s.size());
    for (casadi_int i=0; i<rowperm.size(); ++i) {
      // Permute equations
      daenew[i] = this->dae[rowperm[i]];

      // Permute variables
      snew[i] = this->s[colperm[i]];
      sdotnew[i] = this->sdot[colperm[i]];
    }
    this->dae = daenew;
    this->s = snew;
    this->sdot = sdotnew;

    // Differentiate to write the sorted ODE as a function of state derivatives
    MX J = jacobian(vertcat(this->dae), vertcat(this->sdot));

    // Explicit ODE
    vector<MX> new_ode;

    // Loop over blocks
    for (casadi_int b=0; b<nb; ++b) {

      // Get variables in the block
      vector<MX> xb(this->s.begin()+colblock[b], this->s.begin()+colblock[b+1]);
      vector<MX> xdotb(this->sdot.begin()+colblock[b], this->sdot.begin()+colblock[b+1]);

      // Get equations in the block
      vector<MX> fb(this->dae.begin()+rowblock[b], this->dae.begin()+rowblock[b+1]);

      // Get local Jacobian
      MX Jb = J(Slice(rowblock[b], rowblock[b+1]), Slice(colblock[b], colblock[b+1]));

      // If Jb depends on xb, then the state derivative does not enter linearly
      // in the ODE and we cannot solve for the state derivative
      casadi_assert(!depends_on(Jb, vertcat(xdotb)),
        "Cannot find an explicit expression for variable(s) " + str(xb));

      // Divide fb into a part which depends on vb and a part which doesn't according to
      // "fb == mul(Jb, vb) + fb_res"
      vector<MX> fb_res = substitute(fb, xdotb, vector<MX>(xdotb.size(), 0));

      // Solve for vb
      vector<MX> fb_exp = vertsplit(solve(Jb, -vertcat(fb_res)));

      // Add to explicitly determined equations and variables
      new_ode.insert(new_ode.end(), fb_exp.begin(), fb_exp.end());
    }

    // Eliminate inter-dependencies
    vector<MX> ex;
    substitute_inplace(this->sdot, new_ode, ex, false);

    // Add to explicit differential states and ODE
    this->x.insert(this->x.end(), this->s.begin(), this->s.end());
    this->ode.insert(this->ode.end(), new_ode.begin(), new_ode.end());
    this->dae.clear();
    this->s.clear();
    this->sdot.clear();
  }

  void DaeBuilder::eliminate_alg() {
    // Only works if there are no i
    eliminate_d();

    // Quick return if there are no algebraic states
    if (this->z.empty()) return;

    // Write the algebraic equations as a function of the algebraic states
    Function f("f", {vertcat(this->z)}, {vertcat(this->alg)});

    // Get the sparsity of the Jacobian which can be used to determine which
    // variable can be calculated from which other
    Sparsity sp = f.sparsity_jac(0, 0);
    casadi_assert_dev(sp.is_square());

    // BLT transformation
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    casadi_int nb = sp.btf(rowperm, colperm, rowblock, colblock,
                                  coarse_rowblock, coarse_colblock);

    // Resort equations and variables
    vector<MX> algnew(this->z.size()), znew(this->z.size());
    for (casadi_int i=0; i<rowperm.size(); ++i) {
      // Permute equations
      algnew[i] = this->alg[rowperm[i]];

      // Permute variables
      znew[i] = this->z[colperm[i]];
    }
    this->alg = algnew;
    this->z = znew;

    // Rewrite the sorted algebraic equations as a function of the algebraic states
    f = Function("f", {vertcat(this->z)}, {vertcat(this->alg)});

    // Variables where we have found an explicit expression and where we haven't
    vector<MX> z_exp, z_imp;

    // Explicit and implicit equations
    vector<MX> f_exp, f_imp;

    // Loop over blocks
    for (casadi_int b=0; b<nb; ++b) {

      // Get local variables
      vector<MX> zb(this->z.begin()+colblock[b], this->z.begin()+colblock[b+1]);

      // Get local equations
      vector<MX> fb(this->alg.begin()+rowblock[b], this->alg.begin()+rowblock[b+1]);

      // Get local Jacobian
      MX Jb = jacobian(vertcat(fb), vertcat(zb));

      // If Jb depends on zb, then we cannot (currently) solve for it explicitly
      if (depends_on(Jb, vertcat(zb))) {

        // Add the equations to the new list of algebraic equations
        f_imp.insert(f_imp.end(), fb.begin(), fb.end());

        // ... and the variables accordingly
        z_imp.insert(z_imp.end(), zb.begin(), zb.end());

      } else { // The variables that we wish to determine enter linearly

        // Divide fb into a part which depends on vb and a part which doesn't
        // according to "fb == mul(Jb, vb) + fb_res"
        vector<MX> fb_res = substitute(fb, zb, vector<MX>(zb.size(), 0));

        // Solve for vb
        vector<MX> fb_exp = vertsplit(solve(Jb, -vertcat(fb_res)));

        // Add to explicitly determined equations and variables
        z_exp.insert(z_exp.end(), zb.begin(), zb.end());
        f_exp.insert(f_exp.end(), fb_exp.begin(), fb_exp.end());
      }
    }

    // Eliminate inter-dependencies in fb_exp
    vector<MX> ex;
    substitute_inplace(z_exp, f_exp, ex, false);

    // Add to the beginning of the dependent variables
    // (since the other dependent variable might depend on them)
    this->d.insert(this->d.begin(), z_exp.begin(), z_exp.end());
    this->ddef.insert(this->ddef.begin(), f_exp.begin(), f_exp.end());

    // Save new algebraic equations
    this->z = z_imp;
    this->alg = f_imp;

    // Eliminate new dependent variables from the other equations
    eliminate_d();
  }

  void DaeBuilder::make_explicit() {
    // Only works if there are no i
    eliminate_d();

    // Start by transforming to semi-explicit form
    make_semi_explicit();

    // Then eliminate the algebraic variables
    eliminate_alg();

    // Error if still algebraic variables
    casadi_assert(this->z.empty(), "Failed to eliminate algebraic variables");
  }

  const Variable& DaeBuilder::variable(const std::string& name) const {
    return const_cast<DaeBuilder*>(this)->variable(name);
  }

  Variable& DaeBuilder::variable(const std::string& name) {
    // Find the variable
    VarMap::iterator it = varmap_.find(name);
    if (it==varmap_.end()) {
      casadi_error("No such variable: \"" + name + "\".");
    }

    // Return the variable
    return it->second;
  }

  void DaeBuilder::add_variable(const std::string& name, const Variable& var) {
    // Try to find the component
    if (varmap_.find(name)!=varmap_.end()) {
      stringstream ss;
      casadi_error("Variable \"" + name + "\" has already been added.");
    }

    // Add to the map of all variables
    varmap_[name] = var;
  }

  MX DaeBuilder::add_variable(const std::string& name, casadi_int n) {
    return add_variable(name, Sparsity::dense(n));
  }

  MX DaeBuilder::add_variable(const std::string& name, const Sparsity& sp) {
    Variable v(name, sp);
    add_variable(name, v);
    return v.v;
  }

  MX DaeBuilder::add_x(const std::string& name, casadi_int n) {
    if (name.empty()) return add_x("x" + str(this->x.size()), n);
    MX new_x = add_variable(name, n);
    this->x.push_back(new_x);
    return new_x;
  }

  MX DaeBuilder::add_q(const std::string& name, casadi_int n) {
    if (name.empty()) return add_q("q" + str(this->q.size()), n);
    MX new_q = add_variable(name, n);
    this->q.push_back(new_q);
    return new_q;
  }

  std::pair<MX, MX> DaeBuilder::add_s(const std::string& name, casadi_int n) {
    if (name.empty()) return add_s("s" + str(this->s.size()), n);
    Variable v(name, Sparsity::dense(n));
    add_variable(name, v);
    this->s.push_back(v.v);
    this->sdot.push_back(v.d);
    return std::pair<MX, MX>(v.v, v.d);
  }

  MX DaeBuilder::add_z(const std::string& name, casadi_int n) {
    if (name.empty()) return add_z("z" + str(this->z.size()), n);
    MX new_z = add_variable(name, n);
    this->z.push_back(new_z);
    return new_z;
  }

  MX DaeBuilder::add_p(const std::string& name, casadi_int n) {
    if (name.empty()) return add_p("p" + str(this->p.size()), n);
    MX new_p = add_variable(name, n);
    this->p.push_back(new_p);
    return new_p;
  }

  MX DaeBuilder::add_u(const std::string& name, casadi_int n) {
    if (name.empty()) return add_u("u" + str(this->u.size()), n);
    MX new_u = add_variable(name, n);
    this->u.push_back(new_u);
    return new_u;
  }

  MX DaeBuilder::add_aux(const std::string& name, casadi_int n) {
    if (name.empty()) return add_aux("aux" + str(this->aux.size()), n);
    MX new_aux = add_variable(name, n);
    this->aux.push_back(new_aux);
    return new_aux;
  }

  MX DaeBuilder::add_d(const std::string& name, const MX& new_ddef) {
    MX new_d = add_variable(name, new_ddef.sparsity());
    this->d.push_back(new_d);
    this->ddef.push_back(new_ddef);
    this->lam_ddef.push_back(MX::sym("lam_" + name, new_ddef.sparsity()));
    return new_d;
  }

  MX DaeBuilder::add_y(const std::string& name, const MX& new_ydef) {
    MX new_y = add_variable(name, new_ydef.sparsity());
    this->y.push_back(new_y);
    this->ydef.push_back(new_ydef);
    this->lam_ydef.push_back(MX::sym("lam_" + name, new_ydef.sparsity()));
    return new_y;
  }

  void DaeBuilder::add_ode(const std::string& name, const MX& new_ode) {
    this->ode.push_back(new_ode);
    this->lam_ode.push_back(MX::sym("lam_" + name, new_ode.sparsity()));
  }

  void DaeBuilder::add_dae(const std::string& name, const MX& new_dae) {
    this->dae.push_back(new_dae);
    this->lam_dae.push_back(MX::sym("lam_" + name, new_dae.sparsity()));
  }

  void DaeBuilder::add_alg(const std::string& name, const MX& new_alg) {
    this->alg.push_back(new_alg);
    this->lam_alg.push_back(MX::sym("lam_" + name, new_alg.sparsity()));
  }

  void DaeBuilder::add_quad(const std::string& name, const MX& new_quad) {
    this->quad.push_back(new_quad);
    this->lam_quad.push_back(MX::sym("lam_" + name, new_quad.sparsity()));
  }

  void DaeBuilder::sanity_check() const {
    // Time
    casadi_assert(this->t.is_symbolic(), "Non-symbolic time t");
    casadi_assert(this->t.is_scalar(), "Non-scalar time t");

    // Differential states
    casadi_assert(this->x.size()==this->ode.size(),
                          "x and ode have different lengths");
    for (casadi_int i=0; i<this->x.size(); ++i) {
      casadi_assert(this->x[i].size()==this->ode[i].size(),
                            "ode has wrong dimensions");
      casadi_assert(this->x[i].is_symbolic(), "Non-symbolic state x");
    }

    // DAE
    casadi_assert(this->s.size()==this->sdot.size(),
                          "s and sdot have different lengths");
    casadi_assert(this->s.size()==this->dae.size(),
                          "s and dae have different lengths");
    for (casadi_int i=0; i<this->s.size(); ++i) {
      casadi_assert(this->s[i].is_symbolic(), "Non-symbolic state s");
      casadi_assert(this->s[i].size()==this->sdot[i].size(),
                            "sdot has wrong dimensions");
      casadi_assert(this->s[i].size()==this->dae[i].size(),
                            "dae has wrong dimensions");
    }

    // Algebraic variables/equations
    casadi_assert(this->z.size()==this->alg.size(),
                          "z and alg have different lengths");
    for (casadi_int i=0; i<this->z.size(); ++i) {
      casadi_assert(this->z[i].is_symbolic(), "Non-symbolic algebraic variable z");
      casadi_assert(this->z[i].size()==this->alg[i].size(),
                            "alg has wrong dimensions");
    }

    // Quadrature states/equations
    casadi_assert(this->q.size()==this->quad.size(),
                          "q and quad have different lengths");
    for (casadi_int i=0; i<this->q.size(); ++i) {
      casadi_assert(this->q[i].is_symbolic(), "Non-symbolic quadrature state q");
      casadi_assert(this->q[i].size()==this->quad[i].size(),
                            "quad has wrong dimensions");
    }

    // Intermediate variables
    casadi_assert(this->d.size()==this->ddef.size(),
                          "d and ddef have different lengths");
    for (casadi_int i=0; i<this->d.size(); ++i) {
      casadi_assert(this->d[i].is_symbolic(), "Non-symbolic dependent parameter d");
      casadi_assert(this->d[i].size()==this->ddef[i].size(),
                            "ddef has wrong dimensions");
    }

    // Output equations
    casadi_assert(this->y.size()==this->ydef.size(),
                          "y and ydef have different lengths");
    for (casadi_int i=0; i<this->y.size(); ++i) {
      casadi_assert(this->y[i].is_symbolic(), "Non-symbolic output y");
      casadi_assert(this->y[i].size()==this->ydef[i].size(),
                            "ydef has wrong dimensions");
    }

    // Control
    for (casadi_int i=0; i<this->u.size(); ++i) {
      casadi_assert(this->u[i].is_symbolic(), "Non-symbolic control u");
    }

    // Parameter
    for (casadi_int i=0; i<this->p.size(); ++i) {
      casadi_assert(this->p[i].is_symbolic(), "Non-symbolic parameter p");
    }
  }

  std::string DaeBuilder::qualified_name(const XmlNode& nn) {
    // Stringstream to assemble name
    stringstream qn;

    for (casadi_int i=0; i<nn.size(); ++i) {
      // Add a dot
      if (i!=0) qn << ".";

      // Get the name part
      qn << nn[i].getAttribute("name");

      // Get the index, if any
      if (nn[i].size()>0) {
        casadi_int ind;
        nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText(ind);
        qn << "[" << ind << "]";
      }
    }

    // Return the name
    return qn.str();
  }

  MX DaeBuilder::var(const std::string& name) const {
    return variable(name).v;
  }

  MX DaeBuilder::der(const std::string& name) const {
    return variable(name).d;
  }

  MX DaeBuilder::der(const MX& var) const {
    casadi_assert_dev(var.is_column() && var.is_symbolic());
    MX ret = MX::zeros(var.sparsity());
    for (casadi_int i=0; i<ret.nnz(); ++i) {
      ret.nz(i) = der(var.nz(i).name());
    }
    return ret;
  }

  void DaeBuilder::split_dae() {
    // Only works if there are no d
    eliminate_d();

    // Quick return if no s
    if (this->s.empty()) return;

    // We investigate the interdependencies in sdot -> dae
    Function f("f", {vertcat(this->sdot)}, {vertcat(this->dae)});

    // Number of s
    casadi_int ns = f.nnz_in(0);
    casadi_assert_dev(f.nnz_out(0)==ns);

    // Input/output buffers
    vector<bvec_t> f_sdot(ns, 1);
    vector<bvec_t> f_dae(ns, 0);

    // Propagate to f_dae
    f({get_ptr(f_sdot)}, {get_ptr(f_dae)});

    // Get the new differential and algebraic equations
    vector<MX> new_dae, new_alg;
    for (casadi_int i=0; i<ns; ++i) {
      if (f_dae[i]==bvec_t(1)) {
        new_dae.push_back(this->dae[i]);
      } else {
        casadi_assert_dev(f_dae[i]==bvec_t(0));
        new_alg.push_back(this->dae[i]);
      }
    }

    // Seed all outputs
    std::fill(f_dae.begin(), f_dae.end(), 1);

    // Propagate to f_sdot
    std::fill(f_sdot.begin(), f_sdot.end(), 0);
    f.rev({get_ptr(f_sdot)}, {get_ptr(f_dae)});

    // Get the new algebraic variables and new states
    vector<MX> new_s, new_sdot, new_z;
    for (casadi_int i=0; i<ns; ++i) {
      if (f_sdot[i]==bvec_t(1)) {
        new_s.push_back(this->s[i]);
        new_sdot.push_back(this->sdot[i]);
      } else {
        casadi_assert_dev(f_sdot[i]==bvec_t(0));
        new_z.push_back(this->s[i]);
      }
    }

    // Make sure split was successful
    casadi_assert_dev(new_dae.size()==new_s.size());

    // Divide up the s and dae
    this->dae = new_dae;
    this->s = new_s;
    this->sdot = new_sdot;
    this->alg.insert(this->alg.end(), new_alg.begin(), new_alg.end());
    this->z.insert(this->z.end(), new_z.begin(), new_z.end());
  }

  std::string DaeBuilder::unit(const std::string& name) const {
    return variable(name).unit;
  }

  std::string DaeBuilder::unit(const MX& var) const {
    casadi_assert(!var.is_column() && var.is_valid_input(),
                          "DaeBuilder::unit: Argument must be a symbolic vector");
    if (var.is_empty()) {
      return "n/a";
    } else {
      std::vector<MX> prim = var.primitives();
      string ret = unit(prim.at(0).name());
      for (casadi_int i=1; i<prim.size(); ++i) {
        casadi_assert(ret == unit(prim.at(i).name()),
                              "DaeBuilder::unit: Argument has mixed units");
      }
      return ret;
    }
  }

  void DaeBuilder::set_unit(const std::string& name, const std::string& val) {
    variable(name).unit = val;
  }

  double DaeBuilder::nominal(const std::string& name) const {
    return variable(name).nominal;
  }

  void DaeBuilder::set_nominal(const std::string& name, double val) {
    variable(name).nominal = val;
  }

  std::vector<double> DaeBuilder::nominal(const MX& var) const {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::nominal: Argument must be a symbolic vector");
    std::vector<double> ret(var.nnz());
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      ret[i] = nominal(prim.at(i).name());
    }
    return ret;
  }

  void DaeBuilder::set_nominal(const MX& var, const std::vector<double>& val) {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::nominal: Argument must be a symbolic vector");
    casadi_assert(var.nnz()==var.nnz(), "DaeBuilder::nominal: Dimension mismatch");
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      set_nominal(prim.at(i).name(), val.at(i));
    }
  }

  std::vector<double> DaeBuilder::attribute(getAtt f, const MX& var, bool normalized) const {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::attribute: Argument must be a symbolic vector");
    std::vector<double> ret(var.nnz());
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      ret[i] = (this->*f)(prim[i].name(), normalized);
    }
    return ret;
  }

  MX DaeBuilder::attribute(getAttS f, const MX& var) const {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::attribute: Argument must be a symbolic vector");
    MX ret = MX::zeros(var.sparsity());
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      ret.nz(i) = (this->*f)(prim[i].name());
    }
    return ret;
  }

  void DaeBuilder::set_attribute(setAtt f, const MX& var, const std::vector<double>& val,
                                 bool normalized) {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::set_attribute: Argument must be a symbolic vector");
    casadi_assert(var.nnz()==val.size(), "DaeBuilder::set_attribute: Dimension mismatch");
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      (this->*f)(prim[i].name(), val[i], normalized);
    }
  }

  void DaeBuilder::set_attribute(setAttS f, const MX& var, const MX& val) {
    casadi_assert(var.is_column() && var.is_valid_input(),
                          "DaeBuilder::set_attribute: Argument must be a symbolic vector");
    casadi_assert(var.sparsity()==val.sparsity(),
                          "DaeBuilder::set_attribute: Sparsity mismatch");
    std::vector<MX> prim = var.primitives();
    for (casadi_int i=0; i<prim.size(); ++i) {
      casadi_assert_dev(prim[i].nnz()==1);
      (this->*f)(var.nz(i).name(), val.nz(i));
    }
  }

  double DaeBuilder::min(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.min / v.nominal : v.min;
  }

  std::vector<double> DaeBuilder::min(const MX& var, bool normalized) const {
    return attribute(&DaeBuilder::min, var, normalized);
  }

  void DaeBuilder::set_min(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.min = normalized ? val*v.nominal : val;
  }

  void DaeBuilder::set_min(const MX& var, const std::vector<double>& val, bool normalized) {
    set_attribute(&DaeBuilder::set_min, var, val, normalized);
  }

  double DaeBuilder::max(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.max / v.nominal : v.max;
  }

  std::vector<double> DaeBuilder::max(const MX& var, bool normalized) const {
    return attribute(&DaeBuilder::max, var, normalized);
  }

  void DaeBuilder::set_max(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.max = normalized ? val*v.nominal : val;
  }

  void DaeBuilder::set_max(const MX& var, const std::vector<double>& val, bool normalized) {
    set_attribute(&DaeBuilder::set_max, var, val, normalized);
  }

  double DaeBuilder::guess(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.guess / v.nominal : v.guess;
  }

  std::vector<double> DaeBuilder::guess(const MX& var, bool normalized) const {
    return attribute(&DaeBuilder::guess, var, normalized);
  }

  void DaeBuilder::set_guess(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.guess = normalized ? val*v.nominal : val;
  }

  void DaeBuilder::set_guess(const MX& var, const std::vector<double>& val,
                                    bool normalized) {
    set_attribute(&DaeBuilder::set_guess, var, val, normalized);
  }

  double DaeBuilder::start(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.start / v.nominal : v.start;
  }

  std::vector<double> DaeBuilder::start(const MX& var, bool normalized) const {
    return attribute(&DaeBuilder::start, var, normalized);
  }

  void DaeBuilder::set_start(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.start = normalized ? val*v.nominal : val;
  }

  void DaeBuilder::set_start(const MX& var, const std::vector<double>& val, bool normalized) {
    set_attribute(&DaeBuilder::set_start, var, val, normalized);
  }

  double DaeBuilder::derivative_start(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.derivative_start / v.nominal : v.derivative_start;
  }

  std::vector<double> DaeBuilder::derivative_start(const MX& var, bool normalized) const {
    return attribute(&DaeBuilder::derivative_start, var, normalized);
  }

  void DaeBuilder::set_derivative_start(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.derivative_start = normalized ? val*v.nominal : val;
  }

  void DaeBuilder::set_derivative_start(const MX& var, const std::vector<double>& val,
                                       bool normalized) {
    set_attribute(&DaeBuilder::set_derivative_start, var, val, normalized);
  }

  std::string DaeBuilder::name_in(DaeBuilderIn ind) {
    switch (ind) {
    case DAE_BUILDER_T: return "t";
    case DAE_BUILDER_C: return "c";
    case DAE_BUILDER_P: return "p";
    case DAE_BUILDER_D: return "d";
    case DAE_BUILDER_U: return "u";
    case DAE_BUILDER_X: return "x";
    case DAE_BUILDER_S: return "s";
    case DAE_BUILDER_SDOT: return "sdot";
    case DAE_BUILDER_Z: return "z";
    case DAE_BUILDER_Q: return "q";
    case DAE_BUILDER_W: return "w";
    case DAE_BUILDER_Y: return "y";
    default: return "";
    }
  }

  DaeBuilder::DaeBuilderIn DaeBuilder::enum_in(const std::string& id) {
    if (id=="t") {
      return DAE_BUILDER_T;
    } else if (id=="c") {
      return DAE_BUILDER_C;
    } else if (id=="p") {
      return DAE_BUILDER_P;
    } else if (id=="d") {
      return DAE_BUILDER_D;
    } else if (id=="u") {
      return DAE_BUILDER_U;
    } else if (id=="x") {
      return DAE_BUILDER_X;
    } else if (id=="s") {
      return DAE_BUILDER_S;
    } else if (id=="sdot") {
      return DAE_BUILDER_SDOT;
    } else if (id=="z") {
      return DAE_BUILDER_Z;
    } else if (id=="q") {
      return DAE_BUILDER_Q;
    } else if (id=="w") {
      return DAE_BUILDER_W;
    } else if (id=="y") {
      return DAE_BUILDER_Y;
    } else {
      return DAE_BUILDER_NUM_IN;
    }
  }

  std::vector<DaeBuilder::DaeBuilderIn>
  DaeBuilder::enum_in(const std::vector<std::string>& id) {
    std::vector<DaeBuilderIn> ret(id.size());
    for (casadi_int i=0; i<id.size(); ++i) {
      ret[i] = enum_in(id[i]);
    }
    return ret;
  }

  std::string DaeBuilder::name_out(DaeBuilderOut ind) {
    switch (ind) {
    case DAE_BUILDER_DDEF: return "ddef";
    case DAE_BUILDER_WDEF: return "wdef";
    case DAE_BUILDER_ODE: return "ode";
    case DAE_BUILDER_DAE: return "dae";
    case DAE_BUILDER_ALG: return "alg";
    case DAE_BUILDER_QUAD: return "quad";
    case DAE_BUILDER_YDEF: return "ydef";
    default: return "";
    }
  }

  DaeBuilder::DaeBuilderOut DaeBuilder::enum_out(const std::string& id) {
    if (id=="ddef") {
      return DAE_BUILDER_DDEF;
    } else if (id=="wdef") {
      return DAE_BUILDER_WDEF;
    } else if (id=="ode") {
      return DAE_BUILDER_ODE;
    } else if (id=="dae") {
      return DAE_BUILDER_DAE;
    } else if (id=="alg") {
      return DAE_BUILDER_ALG;
    } else if (id=="quad") {
      return DAE_BUILDER_QUAD;
    } else if (id=="ydef") {
      return DAE_BUILDER_YDEF;
    } else {
      return DAE_BUILDER_NUM_OUT;
    }
  }

  std::vector<DaeBuilder::DaeBuilderOut>
  DaeBuilder::enum_out(const std::vector<std::string>& id) {
    std::vector<DaeBuilderOut> ret(id.size());
    for (casadi_int i=0; i<id.size(); ++i) {
      ret[i] = enum_out(id[i]);
    }
    return ret;
  }

  std::string DaeBuilder::name_in() {
    stringstream ss;
    ss << "[";
    for (casadi_int i=0; i!=DAE_BUILDER_NUM_IN; ++i) {
      if (i!=0) ss << ",";
      ss << name_in(static_cast<DaeBuilderIn>(i));
    }
    ss << "]";
    return ss.str();
  }

  std::string DaeBuilder::name_out() {
    stringstream ss;
    ss << "[";
    for (casadi_int i=0; i!=DAE_BUILDER_NUM_OUT; ++i) {
      if (i!=0) ss << ",";
      ss << name_out(static_cast<DaeBuilderOut>(i));
    }
    ss << "]";
    return ss.str();
  }

  std::vector<MX> DaeBuilder::input(DaeBuilderIn ind) const {
    switch (ind) {
    case DAE_BUILDER_T: return vector<MX>(1, this->t);
    case DAE_BUILDER_C: return this->c;
    case DAE_BUILDER_P: return this->p;
    case DAE_BUILDER_D: return this->d;
    case DAE_BUILDER_U: return this->u;
    case DAE_BUILDER_X: return this->x;
    case DAE_BUILDER_S: return this->s;
    case DAE_BUILDER_SDOT: return this->sdot;
    case DAE_BUILDER_Z: return this->z;
    case DAE_BUILDER_Q: return this->q;
    case DAE_BUILDER_W: return this->w;
    case DAE_BUILDER_Y: return this->y;
    default: return std::vector<MX>();
    }
  }

  std::vector<MX> DaeBuilder::input(std::vector<DaeBuilderIn>& ind) const {
    vector<MX> ret(ind.size());
    for (casadi_int i=0; i<ind.size(); ++i) {
      ret[i] = vertcat(input(ind[i]));
    }
    return ret;
  }

  std::vector<MX> DaeBuilder::output(DaeBuilderOut ind) const {
    switch (ind) {
    case DAE_BUILDER_DDEF: return this->ddef;
    case DAE_BUILDER_WDEF: return this->wdef;
    case DAE_BUILDER_ODE: return this->ode;
    case DAE_BUILDER_DAE: return this->dae;
    case DAE_BUILDER_ALG: return this->alg;
    case DAE_BUILDER_QUAD: return this->quad;
    case DAE_BUILDER_YDEF: return this->ydef;
    default: return std::vector<MX>();
    }
  }

  std::vector<MX> DaeBuilder::output(std::vector<DaeBuilderOut>& ind) const {
    vector<MX> ret(ind.size());
    for (casadi_int i=0; i<ind.size(); ++i) {
      ret[i] = vertcat(output(ind[i]));
    }
    return ret;
  }

  std::vector<MX> DaeBuilder::multiplier(DaeBuilderOut ind) const {
    switch (ind) {
    case DAE_BUILDER_DDEF: return this->lam_ddef;
    case DAE_BUILDER_WDEF: return this->lam_wdef;
    case DAE_BUILDER_ODE: return this->lam_ode;
    case DAE_BUILDER_DAE: return this->lam_dae;
    case DAE_BUILDER_ALG: return this->lam_alg;
    case DAE_BUILDER_QUAD: return this->lam_quad;
    case DAE_BUILDER_YDEF: return this->lam_ydef;
    default: return std::vector<MX>();
    }
  }

  MX DaeBuilder::add_lc(const std::string& name,
                        const std::vector<std::string>& f_out) {
    // Make sure object valid
    sanity_check();

    // Make sure name is valid
    casadi_assert(!name.empty(), "DaeBuilder::add_lc: \"name\" is empty");
    for (string::const_iterator i=name.begin(); i!=name.end(); ++i) {
      casadi_assert(isalnum(*i),
                            "DaeBuilder::add_lc: \"name\" must be alphanumeric");
    }

    // Get a reference to the expression
    MX& ret = lin_comb_[name];
    if (!ret.is_empty()) casadi_warning("DaeBuilder::add_lc: Overwriting " << name);
    ret = 0;

    // Get indices of outputs
    std::vector<DaeBuilderOut> f_out_enum(f_out.size());
    std::vector<bool> in_use(DAE_BUILDER_NUM_OUT, false);
    for (casadi_int i=0; i<f_out.size(); ++i) {
      DaeBuilderOut oind = enum_out(f_out[i]);
      casadi_assert(oind!=DAE_BUILDER_NUM_OUT,
        "DaeBuilder::add_lc: No output expression " + f_out[i] + ". "
        "Valid expressions are " + name_out());
      casadi_assert(!in_use[oind],
        "DaeBuilder::add_lc: Duplicate expression " + f_out[i]);
      in_use[oind] = true;

      // Add linear combination of expressions
      vector<MX> res=output(oind), lam_res=multiplier(oind);
      for (casadi_int i=0; i<res.size(); ++i) {
        ret += dot(lam_res[i], res[i]);
      }
    }

    // Return the (cached) expression
    return ret;
  }

  Function DaeBuilder::create(const std::string& fname,
                                const std::vector<std::string>& s_in,
                                const std::vector<std::string>& s_out) const {
    // Collect function inputs
    vector<MX> ret_in(s_in.size());
    std::vector<bool> input_used(DAE_BUILDER_NUM_IN, false);
    std::vector<bool> output_used(DAE_BUILDER_NUM_IN, false);
    for (vector<string>::const_iterator s_in_it=s_in.begin(); s_in_it!=s_in.end(); ++s_in_it) {
      // Primal variable
      DaeBuilderIn iind = enum_in(*s_in_it);
      if (iind!=DAE_BUILDER_NUM_IN) {
        casadi_assert(!input_used[iind],
          "DaeBuilder::function: Duplicate expression " + *s_in_it);
        input_used[iind] = true;
        ret_in[s_in_it-s_in.begin()] = vertcat(input(iind));
        continue;
      }

      // Dual variable
      if (s_in_it->size()>4 && s_in_it->substr(0, 4)=="lam_") {
        DaeBuilderOut oind = enum_out(s_in_it->substr(4, string::npos));
        if (oind!=DAE_BUILDER_NUM_OUT) {
          casadi_assert(!output_used[oind],
            "DaeBuilder::function: Duplicate expression " + *s_in_it);
          output_used[oind] = true;
          ret_in[s_in_it-s_in.begin()] = vertcat(multiplier(oind));
          continue;
        }
      }

      // Error if reached this point
      stringstream ss;
      ss << "DaeBuilder::function: No input expression " << *s_in_it << "." << endl;
      ss << "Valid expressions are: [";
      for (casadi_int i=0; i!=DAE_BUILDER_NUM_IN; ++i) {
        if (i!=0) ss << ", ";
        ss << name_in(static_cast<DaeBuilderIn>(i));
      }
      for (casadi_int i=0; i!=DAE_BUILDER_NUM_OUT; ++i) {
        ss << ", lam_" << name_out(static_cast<DaeBuilderOut>(i));
      }
      ss << "]";
      casadi_error(ss.str());
    }

    // Function outputs
    vector<MX> ret_out(s_out.size());
    vector<bool> assigned(s_out.size(), false);

    // List of valid attributes
    enum Attributes {
      ATTR_TRANSPOSE,
      ATTR_TRIU,
      ATTR_TRIL,
      ATTR_DENSIFY};

    // Separarate attributes
    vector<vector<Attributes> > attr(s_out.size());
    vector<string> s_out_noatt = s_out;
    for (casadi_int i=0; i<s_out_noatt.size(); ++i) {
      // Currently processed string
      string& s = s_out_noatt[i];

      // Loop over attributes
      while (true) {
        // Find the first underscore separator
        size_t pos = s.find('_');
        if (pos>=s.size()) break; // No more underscore
        string a = s.substr(0, pos);

        // Abort if not attribute
        if (a=="transpose") {
          attr[i].push_back(ATTR_TRANSPOSE);
        } else if (a=="triu") {
          attr[i].push_back(ATTR_TRIU);
        } else if (a=="tril") {
          attr[i].push_back(ATTR_TRIL);
        } else if (a=="densify") {
          attr[i].push_back(ATTR_DENSIFY);
        } else {
          // No more attribute
          break;
        }

        // Strip attribute
        s = s.substr(pos+1, string::npos);
      }
    }

    // Non-differentiated outputs
    fill(output_used.begin(), output_used.end(), false);
    for (casadi_int i=0; i<s_out_noatt.size(); ++i) {
      DaeBuilderOut oind = enum_out(s_out_noatt[i]);
      if (oind!=DAE_BUILDER_NUM_OUT) {
        casadi_assert(!output_used[oind],
          "DaeBuilder::function: Duplicate expression " + s_out_noatt[i]);
        output_used[oind] = true;
        ret_out[i] = vertcat(output(oind));
        assigned[i] = true;
      }
    }

    // Linear combination of outputs
    for (casadi_int i=0; i<s_out_noatt.size(); ++i) {
      if (assigned[i]) continue;
      std::map<std::string, MX>::const_iterator j=lin_comb_.find(s_out_noatt[i]);
      if (j!=lin_comb_.end()) {
        ret_out[i] = j->second;
        assigned[i] = true;
      }
    }

    // Determine which Jacobian blocks to generate
    vector<vector<casadi_int> > wanted(DAE_BUILDER_NUM_OUT,
      vector<casadi_int>(DAE_BUILDER_NUM_IN, -1));
    for (casadi_int i=0; i<s_out_noatt.size(); ++i) {
      if (assigned[i]) continue;

      // Get the string without attributes
      const string& so = s_out_noatt[i];

      // Find the first underscore separator
      size_t pos = so.find('_');
      if (pos>=so.size()) continue;

      // Get operation
      string s = so.substr(0, pos);
      if (s!="jac") continue;

      // Get expression to be differentiated
      size_t pos1 = so.find('_', pos+1);
      if (pos1>=so.size()) continue;
      s = so.substr(pos+1, pos1-pos-1);
      DaeBuilderOut oind = enum_out(s);
      if (oind==DAE_BUILDER_NUM_OUT) continue;

      // Jacobian with respect to what variable
      s = so.substr(pos1+1, string::npos);
      DaeBuilderIn iind = enum_in(s);
      if (iind==DAE_BUILDER_NUM_IN) continue;

      // Check if duplicate
      casadi_assert(wanted[oind][iind]==-1,
        "DaeBuilder::function: Duplicate Jacobian " + so);
      wanted[oind][iind] = i;
    }

    // Generate Jacobian blocks
    for (casadi_int oind=0; oind!=DAE_BUILDER_NUM_OUT; ++oind) {
      for (casadi_int iind=0; iind!=DAE_BUILDER_NUM_IN; ++iind) {
        // Skip if not wanted
        if (wanted[oind][iind]==-1) continue;

        // List of blocks to be calculated together, starting with current
        vector<DaeBuilderIn> ib(1, static_cast<DaeBuilderIn>(iind));
        vector<DaeBuilderOut> ob(1, static_cast<DaeBuilderOut>(oind));

        // Add other blocks that can be calculated with the same inputs
        // (typically cheap if forward mode used)
        for (casadi_int oind1=oind+1; oind1!=DAE_BUILDER_NUM_OUT; ++oind1) {
          if (wanted[oind1][iind]>=0) {
            ob.push_back(static_cast<DaeBuilderOut>(oind1));
          }
        }

        // Add other blocks that can be calculated with the same outputs
        // (typically cheap if reverse mode used)
        for (casadi_int iind1=iind+1; iind1!=DAE_BUILDER_NUM_IN; ++iind1) {
          // Do we really want _all_ the input/output combinations?
          bool all_wanted = true;
          for (casadi_int k=0; k<ob.size() && all_wanted; ++k) {
            all_wanted = wanted[ob[k]][iind1]>=0;
          }

          // Add block(s)
          if (all_wanted) {
            ib.push_back(static_cast<DaeBuilderIn>(iind1));
          }
        }

        // When we know which blocks we are interested in, we form the Jacobian
        vector<MX> arg=input(ib), res=output(ob);
        MX J = jacobian(vertcat(res), vertcat(arg));

        // Divide into blocks and copy to output
        vector<vector<MX> > J_all = blocksplit(J, offset(res), offset(arg));
        for (casadi_int ki=0; ki<ib.size(); ++ki) {
          for (casadi_int ko=0; ko<ob.size(); ++ko) {
            casadi_int& ind=wanted[ob[ko]][ib[ki]];
            ret_out[ind] = J_all[ko][ki];
            assigned[ind] = true;
            ind = -1;
          }
        }
      }
    }

    // For all linear combinations
    for (std::map<std::string, MX>::const_iterator lin_comb_it=lin_comb_.begin();
         lin_comb_it!=lin_comb_.end(); ++lin_comb_it) {
      // Determine which Hessian blocks to generate
      wanted.resize(DAE_BUILDER_NUM_IN);
      fill(wanted.begin(), wanted.end(), vector<casadi_int>(DAE_BUILDER_NUM_IN, -1));

      for (casadi_int i=0; i<s_out_noatt.size(); ++i) {
        if (assigned[i]) continue;

        // Get the string without attributes
        const string& so = s_out_noatt[i];

        // Get operation
        size_t pos = so.find('_');
        if (pos>=so.size()) continue;
        string s = so.substr(0, pos);
        if (s!="hes") continue;

        // Get expression to be differentiated
        size_t pos1 = so.find('_', pos+1);
        if (pos1>=so.size()) continue;
        s = so.substr(pos+1, pos1-pos-1);
        if (s!=lin_comb_it->first) continue;

        // Get first derivative
        pos = so.find('_', pos1+1);
        if (pos>=so.size()) continue;
        s = so.substr(pos1+1, pos-pos1-1);
        DaeBuilderIn iind1 = enum_in(s);
        if (iind1==DAE_BUILDER_NUM_IN) continue;

        // Get second derivative
        s = so.substr(pos+1, string::npos);
        DaeBuilderIn iind2 = enum_in(s);
        if (iind2==DAE_BUILDER_NUM_IN) continue;

        // Check if duplicate
        casadi_assert(wanted[iind1][iind2]==-1,
          "DaeBuilder::function: Duplicate Hessian " + so);
        wanted[iind1][iind2] = i;
      }

      // Created wanted Hessian blocks
      for (casadi_int iind1=0; iind1!=DAE_BUILDER_NUM_IN; ++iind1) {
        for (casadi_int iind2=0; iind2!=DAE_BUILDER_NUM_IN; ++iind2) {
          // Skip if not wanted
          if (wanted[iind1][iind2]==-1) continue;

          // List of blocks to be calculated together, starting with current
          vector<DaeBuilderIn> ib1(1, static_cast<DaeBuilderIn>(iind1));
          vector<DaeBuilderIn> ib2(1, static_cast<DaeBuilderIn>(iind2));

          // Add other blocks vertically
          for (casadi_int iind=iind1+1; iind!=DAE_BUILDER_NUM_IN; ++iind) {
            if (wanted[iind][iind2]>=0 || wanted[iind2][iind]>=0) {
              ib1.push_back(static_cast<DaeBuilderIn>(iind));
            }
          }

          // Add other blocks horizontally
          for (casadi_int iind=iind2+1; iind!=DAE_BUILDER_NUM_IN; ++iind) {
            // Do we really want _all_ the blocks?
            bool all_wanted = true;
            for (casadi_int k=0; k<ib1.size() && all_wanted; ++k) {
              all_wanted = wanted[ib1[k]][iind]>=0 || wanted[iind][ib1[k]]>=0;
            }

            // Add block(s)
            if (all_wanted) {
              ib2.push_back(static_cast<DaeBuilderIn>(iind));
            }
          }

          // Symmetric or not?
          bool symmetric=ib1.size()==ib2.size();
          for (casadi_int i=0; symmetric && i<ib1.size(); ++i) symmetric=ib1[i]==ib2[i];

          // Calculate blocks
          vector<vector<MX> > H_all;
          if (symmetric) {
            vector<MX> arg=input(ib1);
            MX H = hessian(lin_comb_it->second, vertcat(arg));
            H_all = blocksplit(H, offset(arg), offset(arg));
          } else {
            vector<MX> arg1=input(ib1), arg2=input(ib2);
            MX g = gradient(lin_comb_it->second, vertcat(arg1));
            MX J = jacobian(g, vertcat(arg2));
            H_all = blocksplit(J, offset(arg1), offset(arg2));
          }
          // Fetch the requested blocks
          for (casadi_int k1=0; k1<ib1.size(); ++k1) {
            for (casadi_int k2=0; k2<ib2.size(); ++k2) {
              casadi_int& ind=wanted[ib1[k1]][ib2[k2]];
              if (ind>=0) {
                ret_out[ind] = H_all[k1][k2];
                assigned[ind] = true;
                ind = -1;
              }
            }
          }
        }
      }
    }

    // Check and post-process outputs
    for (casadi_int i=0; i<s_out.size(); ++i) {
      // Make sure all outputs have been assigned
      if (!assigned[i]) {
        casadi_error("DaeBuilder::function: Cannot treat output expression " + s_out[i]);
      }

      // Apply attributes starting from the right-most one
      MX& r = ret_out[i];
      for (auto a=attr[i].rbegin(); a!=attr[i].rend(); ++a) {
        switch (*a) {
        case ATTR_TRANSPOSE: r = r.T(); break;
        case ATTR_TRIU: r = triu(r); break;
        case ATTR_TRIL: r = tril(r); break;
        case ATTR_DENSIFY: r = densify(r); break;
        }
      }
    }

    // Generate the constructed function
    return Function(fname, ret_in, ret_out, s_in, s_out);
  }

  Function DaeBuilder::add_fun(const Function& f) {
    casadi_assert(!has_fun(f.name()), "Function '" + f.name() + "' already exists");
    fun_.push_back(f);
    return f;
  }

  Function DaeBuilder::add_fun(const std::string& name,
                               const std::vector<std::string>& arg,
                               const std::vector<std::string>& res,
                               const Dict& opts) {
    casadi_assert(!has_fun(name), "Function '" + name + "' already exists");

    // Get inputs
    vector<MX> arg_ex, res_ex;
    for (auto&& s : arg) arg_ex.push_back(var(s));
    for (auto&& s : res) {
      // Find the binding expression FIXME(@jaeandersson)
      casadi_int d_ind;
      for (d_ind=0; d_ind<this->d.size(); ++d_ind) {
        if (s==this->d.at(d_ind).name()) {
          res_ex.push_back(this->ddef.at(d_ind));
          break;
        }
      }
      casadi_assert(d_ind<this->d.size(), "Cannot find dependent '" + s + "'");
    }
    Function ret(name, arg_ex, res_ex, arg, res, opts);
    return add_fun(ret);
  }

  Function DaeBuilder::add_fun(const std::string& name, const Importer& compiler,
                               const Dict& opts) {
    casadi_assert(!has_fun(name), "Function '" + name + "' already exists");
    return add_fun(external(name, compiler, opts));
  }

  bool DaeBuilder::has_fun(const std::string& name) const {
    for (const Function& f : fun_) {
      if (f.name()==name) return true;
    }
    return false;
  }

  Function DaeBuilder::fun(const std::string& name) const {
    casadi_assert(has_fun(name), "No such function: '" + name + "'");
    for (const Function& f : fun_) {
      if (f.name()==name) return f;
    }
    return Function();
  }

} // namespace casadi
