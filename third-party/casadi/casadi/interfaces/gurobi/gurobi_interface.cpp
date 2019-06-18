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


#include "gurobi_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_GUROBI_EXPORT
  casadi_register_conic_gurobi(Conic::Plugin* plugin) {
    plugin->creator = GurobiInterface::creator;
    plugin->name = "gurobi";
    plugin->doc = GurobiInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &GurobiInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_GUROBI_EXPORT casadi_load_conic_gurobi() {
    Conic::registerPlugin(casadi_register_conic_gurobi);
  }

  GurobiInterface::GurobiInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  GurobiInterface::~GurobiInterface() {
    clear_mem();
  }

  Options GurobiInterface::options_
  = {{&Conic::options_},
     {{"vtype",
       {OT_STRINGVECTOR,
        "Type of variables: [CONTINUOUS|binary|integer|semicont|semiint]"}},
      {"gurobi",
       {OT_DICT,
        "Options to be passed to gurobi."}}
     }
  };

  void GurobiInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    std::vector<std::string> vtype;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="vtype") {
        vtype = op.second;
      } else if (op.first=="gurobi") {
        opts_ = op.second;
      }
    }

    // Variable types
    if (!vtype.empty()) {
      casadi_assert(vtype.size()==nx_, "Option 'vtype' has wrong length");
      vtype_.resize(nx_);
      for (casadi_int i=0; i<nx_; ++i) {
        if (vtype[i]=="continuous") {
          vtype_[i] = GRB_CONTINUOUS;
        } else if (vtype[i]=="binary") {
          vtype_[i] = GRB_BINARY;
        } else if (vtype[i]=="integer") {
          vtype_[i] = GRB_INTEGER;
        } else if (vtype[i]=="semicont") {
          vtype_[i] = GRB_SEMICONT;
        } else if (vtype[i]=="semiint") {
          vtype_[i] = GRB_SEMIINT;
        } else {
          casadi_error("No such variable type: " + vtype[i]);
        }
      }
    }

    Sparsity qsum = reshape(sum2(Q_), np_, np_);

    // Block detection
    Sparsity aggregate = qsum+P_;

    std::vector<casadi_int> p;
    casadi_int nb = aggregate.scc(p, r_);

    std::string pattern_message = "Pattern not recognised";

    casadi_assert(p==range(p.size()), pattern_message);

    const casadi_int* row = aggregate.row();
    const casadi_int* colind = aggregate.colind();

    // Check fishbone-structure
    for (casadi_int i=0;i<nb;++i) {
      casadi_int block_size = r_[i+1]-r_[i];
      // number of nonzeros in column r_[i+1]-1
      casadi_int nz = colind[r_[i+1]]-colind[r_[i+1]-1];
      // Last column of block should be dense
      casadi_assert(nz==block_size, pattern_message);
      for (casadi_int k=0;k<block_size-1;++k) {
        casadi_assert(colind[r_[i]+k+1]-colind[r_[i]+k], pattern_message);
        casadi_assert(*(row++)==k+r_[i], pattern_message);
        casadi_assert(*(row++)==r_[i]+block_size-1, pattern_message);
      }

      for (casadi_int k=0;k<block_size;++k)
        casadi_assert(*(row++)==k+r_[i], pattern_message);
    }

    /**
      general soc constraints:
      ||Ax + b ||_2 <= c'x + d

      Need to represent as
      X'X <= Z'Z

      with X and Z helper variables and
      Ax  + b = X
      c'x + d = Z

      [A;c'] x - [X;Z] = [b;d]

      we look for the vertical concatenation of these constraints for all blocks:

      Q(map_Q) [x;X1;Z1;X2;Z2;...] = P.nz(map_P)

    */

    /*

    Aggregate pattern:

    x (x)
     x(x)
    xx(x)
       x (x)
        x(x)
       xx(x)

    We are interested in the parts in parenthesis (target).

    Find out which rows in Q correspond to those targets
    */

    // Lookup vector for target start and end
    std::vector<casadi_int> target_start(nb), target_stop(nb);
    for (casadi_int i=0;i<nb;++i) {
      target_start[i] = (r_[i+1]-1)*aggregate.size1()+r_[i];
      target_stop[i] = target_start[i]+r_[i+1]-r_[i];
    }

    // Collect the nonzero indices in Q that that correspond to the target area
    std::vector<casadi_int> q_nz;
    // Triplet form for map_Q sparsity
    std::vector<casadi_int> q_row, q_col;

    // Loop over Q's columns (decision variables)
    for (casadi_int j=0; j<Q_.size2(); ++j) {
      casadi_int block_index = 0;
      // Loop over Q's rows
      for (casadi_int k=Q_.colind(j); k<Q_.colind(j+1); ++k) {
        casadi_int i = Q_.row(k);

        // Increment block_index if i runs ahead
        while (i>target_stop[block_index] && block_index<nb-1) block_index++;

        if (i>=target_start[block_index] && i<target_stop[block_index]) {
          // Got a nonzero in the target region
          q_nz.push_back(k);
          q_row.push_back(r_[block_index]+i-target_start[block_index]);
          q_col.push_back(j);
        }
      }
    }

    map_Q_ = IM::triplet(q_row, q_col, q_nz, r_[nb], nx_);

    // Add the [X1;Z1;X2;Z2;...] part
    map_Q_ = horzcat(map_Q_, -IM::eye(r_[nb])).T();

    // Get maximum nonzero count of any column
    casadi_int max_nnz = 0;
    for (casadi_int i=0;i<map_Q_.size2();++i) {
      max_nnz = std::max(max_nnz, map_Q_.colind(i+1)-map_Q_.colind(i));
    }

    // ind/val size needs to cover max nonzero count
    indval_size_ = std::max(nx_, max_nnz);

    // Collect the indices for the P target area
    map_P_.resize(r_[nb], -1);
    for (casadi_int i=0;i<nb;++i) {
      for (casadi_int k=P_.colind(r_[i+1]-1); k<P_.colind(r_[i+1]); ++k) {
        casadi_int r = P_.row(k);
        map_P_[r] = k;
      }
    }

    // ind/val size needs to cover blocksize
    for (casadi_int i=0;i<nb;++i)
      indval_size_ = std::max(indval_size_, r_[i+1]-r_[i]);

    // Get the transpose and mapping
    AT_ = A_.transpose(A_mapping_);

    // Temporary memory
    alloc_w(indval_size_, true); // val
    alloc_iw(indval_size_, true); // ind
    alloc_iw(nx_, true); // ind2
  }

  int GurobiInterface::init_mem(void* mem) const {
    auto m = static_cast<GurobiMemory*>(mem);

    // Load environment
    casadi_int flag = GRBloadenv(&m->env, nullptr); // no log file
    casadi_assert(!flag && m->env, "Failed to create GUROBI environment. Flag: "+ str(flag));

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    switch (status) {
    case GRB_LOADED:
      return "LOADED";
    case GRB_OPTIMAL:
      return "OPTIMAL";
    case GRB_INFEASIBLE:
      return "INFEASIBLE";
    case GRB_INF_OR_UNBD:
      return "INF_OR_UNBD";
    case GRB_UNBOUNDED:
      return "UNBOUNDED";
    case GRB_CUTOFF:
      return "CUTOFF";
    case GRB_ITERATION_LIMIT:
      return "ITERATION_LIMIT";
    case GRB_NODE_LIMIT:
      return "NODE_LIMIT";
    case GRB_TIME_LIMIT:
      return "TIME_LIMIT";
    case GRB_SOLUTION_LIMIT:
      return "SOLUTION_LIMIT";
    case GRB_INTERRUPTED:
      return "INTERRUPTED";
    case GRB_NUMERIC:
      return "NUMERIC";
    case GRB_SUBOPTIMAL:
      return "SUBOPTIMAL";
    case GRB_INPROGRESS:
      return "INPROGRESS";
    }
    return "Unknown";
  }

  int GurobiInterface::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<GurobiMemory*>(mem);

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();

    m->fstats.at("preprocessing").tic();

    // Problem has not been solved at this point
    m->success = false;
    m->return_status = -1;

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Inputs
    const double *h=arg[CONIC_H],
      *g=arg[CONIC_G],
      *a=arg[CONIC_A],
      *lba=arg[CONIC_LBA],
      *uba=arg[CONIC_UBA],
      *lbx=arg[CONIC_LBX],
      *ubx=arg[CONIC_UBX],
      *p=arg[CONIC_P],
      *q=arg[CONIC_Q];
      //*x0=arg[CONIC_X0],
      //*lam_x0=arg[CONIC_LAM_X0];

    // Outputs
    double *x=res[CONIC_X],
      *cost=res[CONIC_COST],
      *lam_a=res[CONIC_LAM_A],
      *lam_x=res[CONIC_LAM_X];

    // Temporary memory
    double *val=w; w+=indval_size_;
    int *ind=reinterpret_cast<int*>(iw); iw+=indval_size_;
    int *ind2=reinterpret_cast<int*>(iw); iw+=nx_;

    // Greate an empty model
    GRBmodel *model = nullptr;
    try {
      casadi_int flag = GRBnewmodel(m->env, &model, name_.c_str(), 0,
        nullptr, nullptr, nullptr, nullptr, nullptr);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      // Add variables
      for (casadi_int i=0; i<nx_; ++i) {
        // Get bounds
        double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
        if (isinf(lb)) lb = -GRB_INFINITY;
        if (isinf(ub)) ub =  GRB_INFINITY;

        // Get variable type
        char vtype;
        if (!vtype_.empty()) {
          // Explicitly set 'vtype' takes precedence
          vtype = vtype_.at(i);
        } else if (!discrete_.empty() && discrete_.at(i)) {
          // Variable marked as discrete (integer or binary)
          vtype = lb==0 && ub==1 ? GRB_BINARY : GRB_INTEGER;
        } else {
          // Continious variable
          vtype = GRB_CONTINUOUS;
        }

        // Pass to model
        flag = GRBaddvar(model, 0, nullptr, nullptr, g ? g[i] : 0., lb, ub, vtype, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // Add helper variables for SOCP
      for (casadi_int i=0;i<r_.size()-1;++i) {
        for (casadi_int k=0;k<r_[i+1]-r_[i]-1;++k) {
          flag = GRBaddvar(model, 0, nullptr, nullptr, 0, -GRB_INFINITY, GRB_INFINITY,
                           GRB_CONTINUOUS, nullptr);
          casadi_assert(!flag, GRBgeterrormsg(m->env));
        }
        flag = GRBaddvar(model, 0, nullptr, nullptr, 0, 0, GRB_INFINITY, GRB_CONTINUOUS, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      flag = GRBupdatemodel(model);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      // Add quadratic terms
      const casadi_int *H_colind=H_.colind(), *H_row=H_.row();
      for (int i=0; i<nx_; ++i) {

        // Quadratic term nonzero indices
        casadi_int numqnz = H_colind[1]-H_colind[0];
        for (casadi_int k=0;k<numqnz;++k) ind[k]=H_row[k];
        H_colind++;
        H_row += numqnz;

        // Corresponding column
        casadi_fill(ind2, numqnz, i);

        // Quadratic term nonzeros
        if (h) {
          casadi_copy(h, numqnz, val);
          casadi_scal(numqnz, 0.5, val);
          h += numqnz;
        } else {
          casadi_fill(val, numqnz, 0.);
        }

        // Pass to model
        flag = GRBaddqpterms(model, numqnz, ind, ind2, val);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // Add constraints
      const casadi_int *AT_colind=AT_.colind(), *AT_row=AT_.row();
      for (casadi_int i=0; i<na_; ++i) {
        // Get bounds
        double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;

        casadi_int numnz = 0;
        // Loop over rows
        for (casadi_int k=AT_colind[i]; k<AT_colind[i+1]; ++k) {
          casadi_int j = AT_row[k];

          ind[numnz] = j;
          val[numnz] = a ? a[A_mapping_[k]]  : 0;

          numnz++;
        }

        // Pass to model
        if (isinf(lb)) {
          if (isinf(ub)) {
            // Neither upper or lower bounds, skip
          } else {
            // Only upper bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_LESS_EQUAL, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          }
        } else {
          if (isinf(ub)) {
            // Only lower bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_GREATER_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          } else if (lb==ub) {
            // Upper and lower bounds equal
            flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          } else {
            // Both upper and lower bounds
            flag = GRBaddrangeconstr(model, numnz, ind, val, lb, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          }
        }
      }

      // SOCP helper constraints
      const Sparsity& sp = map_Q_.sparsity();
      const casadi_int* colind = sp.colind();
      const casadi_int* row = sp.row();
      const casadi_int* data = map_Q_.ptr();

      // Loop over columns
      for (casadi_int i=0; i<sp.size2(); ++i) {

        casadi_int numnz = 0;
        // Loop over rows
        for (casadi_int k=colind[i]; k<colind[i+1]; ++k) {
          casadi_int j = row[k];

          ind[numnz] = j;
          val[numnz] = (q && j<nx_) ? q[data[k]] : -1;

          numnz++;
        }

        // Get bound
        double bound = map_P_[i]==-1 ? 0 : -p[map_P_[i]];

        flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, bound, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // Loop over blocks
      for (casadi_int i=0; i<r_.size()-1; ++i) {
        casadi_int block_size = r_[i+1]-r_[i];

        // Indicate x'x - y^2 <= 0
        for (casadi_int j=0;j<block_size;++j) {
          ind[j] = nx_ + r_[i] + j;
          val[j] = j<block_size-1 ? 1 : -1;
        }

        flag = GRBaddqconstr(model, 0, nullptr, nullptr,
          block_size, ind, ind, val,
          GRB_LESS_EQUAL, 0, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      flag = 0;
      for (auto && op : opts_) {
        int ret = GRBgetparamtype(m->env, op.first.c_str());
        switch (ret) {
          case -1:
            casadi_error("Parameter '" + op.first + "' unknown to Gurobi.");
          case 1:
            {
              flag = GRBsetintparam(GRBgetenv(model), op.first.c_str(), op.second);
              break;
            }
          case 2:
              flag = GRBsetdblparam(GRBgetenv(model), op.first.c_str(), op.second);
              break;
          case 3:
            {
              std::string s = op.second;
              flag = GRBsetstrparam(GRBgetenv(model), op.first.c_str(), s.c_str());
              break;
            }
          default:
            casadi_error("Not implememented : " + str(ret));
        }
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();

      // Solve the optimization problem
      flag = GRBoptimize(model);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      m->fstats.at("solver").toc();
      m->fstats.at("postprocessing").tic();

      int optimstatus;
      flag = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      if (verbose_) uout() << "return status: " << return_status_string(optimstatus) <<
        " (" << optimstatus <<")" << std::endl;

      m->return_status = optimstatus;
      m->success = optimstatus==GRB_OPTIMAL;

      // Get the objective value, if requested
      if (cost) {
        flag = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, cost);
        if (flag) cost[0] = casadi::nan;
      }

      // Get the optimal solution, if requested
      if (x) {
        flag = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nx_, x);
        if (flag) fill_n(x, nx_, casadi::nan);
      }

      if (lam_x) fill_n(lam_x, nx_, casadi::nan);
      if (lam_a) fill_n(lam_a, na_, casadi::nan);

      // Free memory
      GRBfreemodel(model);
      m->fstats.at("postprocessing").toc();

    } catch (...) {
      // Free memory
      if (model) GRBfreemodel(model);
      throw;
    }

    // Show statistics
    if (print_time_)  print_fstats(static_cast<ConicMemory*>(mem));
    return 0;
  }

  Dict GurobiInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<GurobiMemory*>(mem);
    stats["return_status"] = return_status_string(m->return_status);
    stats["success"] = m->success;
    return stats;
  }

  GurobiMemory::GurobiMemory() {
    this->env = nullptr;
  }

  GurobiMemory::~GurobiMemory() {
    if (this->env) GRBfreeenv(this->env);
  }

} // namespace casadi
