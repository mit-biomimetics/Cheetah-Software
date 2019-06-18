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

#include "mx_function.hpp"
#include "casadi_misc.hpp"
#include "casadi_common.hpp"
#include "global_options.hpp"
#include "casadi_interrupt.hpp"
#include "io_instruction.hpp"

#include <stack>
#include <typeinfo>

// Throw informative error message
#define CASADI_THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in MXFunction::" FNAME " at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

using namespace std;

namespace casadi {

  MXFunction::MXFunction(const std::string& name,
                         const std::vector<MX>& inputv,
                         const std::vector<MX>& outputv,
                         const std::vector<std::string>& name_in,
                         const std::vector<std::string>& name_out) :
    XFunction<MXFunction, MX, MXNode>(name, inputv, outputv, name_in, name_out) {
  }

  MXFunction::~MXFunction() {
  }

  Options MXFunction::options_
  = {{&FunctionInternal::options_},
     {{"default_in",
       {OT_DOUBLEVECTOR,
        "Default input values"}},
      {"live_variables",
       {OT_BOOL,
        "Reuse variables in the work vector"}}
     }
  };

  MX MXFunction::instruction_MX(casadi_int k) const {
    return algorithm_.at(k).data;
  }

  std::vector<casadi_int> MXFunction::instruction_input(casadi_int k) const {
    auto e = algorithm_.at(k);
    if (e.op==OP_INPUT) {
      const IOInstruction* io = static_cast<const IOInstruction*>(e.data.get());
      return { io->ind() };
    } else {
      return e.arg;
    }
  }

  std::vector<casadi_int> MXFunction::instruction_output(casadi_int k) const {
    auto e = algorithm_.at(k);
    if (e.op==OP_OUTPUT) {
      const IOInstruction* io = static_cast<const IOInstruction*>(e.data.get());
      return { io->ind() };
    } else {
      return e.res;
    }
  }

  void MXFunction::init(const Dict& opts) {
    // Call the init function of the base class
    XFunction<MXFunction, MX, MXNode>::init(opts);
    if (verbose_) casadi_message(name_ + "::init");

    // Default (temporary) options
    bool live_variables = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="default_in") {
        default_in_ = op.second;
      } else if (op.first=="live_variables") {
        live_variables = op.second;
      }
    }

    // Check/set default inputs
    if (default_in_.empty()) {
      default_in_.resize(n_in_, 0);
    } else {
      casadi_assert(default_in_.size()==n_in_,
                            "Option 'default_in' has incorrect length");
    }

    // Stack used to sort the computational graph
    stack<MXNode*> s;

    // All nodes
    vector<MXNode*> nodes;

    // Add the list of nodes
    for (casadi_int ind=0; ind<out_.size(); ++ind) {
      // Loop over primitives of each output
      vector<MX> prim = out_[ind].primitives();
      casadi_int nz_offset=0;
      for (casadi_int p=0; p<prim.size(); ++p) {
        // Get the nodes using a depth first search
        s.push(prim[p].get());
        sort_depth_first(s, nodes);
        // Add an output instruction ("data" below will take ownership)
        nodes.push_back(new Output(prim[p], ind, p, nz_offset));
        // Update offset
        nz_offset += prim[p].nnz();
      }
    }

    // Set the temporary variables to be the corresponding place in the sorted graph
    for (casadi_int i=0; i<nodes.size(); ++i) {
      nodes[i]->temp = i;
    }

    // Place in the algorithm for each node
    vector<casadi_int> place_in_alg;
    place_in_alg.reserve(nodes.size());

    // Input instructions
    vector<pair<casadi_int, MXNode*> > symb_loc;

    // Count the number of times each node is used
    vector<casadi_int> refcount(nodes.size(), 0);

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size());
    for (MXNode* n : nodes) {

      // Get the operation
      casadi_int op = n->op();

      // Store location if parameter (or input)
      if (op==OP_PARAMETER) {
        symb_loc.push_back(make_pair(algorithm_.size(), n));
      }

      // If a new element in the algorithm needs to be added
      if (op>=0) {
        AlgEl ae;
        ae.op = op;
        ae.data.own(n);
        ae.arg.resize(n->n_dep());
        for (casadi_int i=0; i<n->n_dep(); ++i) {
          ae.arg[i] = n->dep(i)->temp;
        }
        ae.res.resize(n->nout());
        if (n->has_output()) {
          fill(ae.res.begin(), ae.res.end(), -1);
        } else if (!ae.res.empty()) {
          ae.res[0] = n->temp;
        }

        // Increase the reference count of the dependencies
        for (casadi_int c=0; c<ae.arg.size(); ++c) {
          if (ae.arg[c]>=0) {
            refcount[ae.arg[c]]++;
          }
        }

        // Save to algorithm
        place_in_alg.push_back(algorithm_.size());
        algorithm_.push_back(ae);

      } else { // Function output node
        // Get the output index
        casadi_int oind = n->which_output();

        // Get the index of the parent node
        casadi_int pind = place_in_alg[n->dep(0)->temp];

        // Save location in the algorithm element corresponding to the parent node
        casadi_int& otmp = algorithm_[pind].res.at(oind);
        if (otmp<0) {
          otmp = n->temp; // First time this function output is encountered, save to algorithm
        } else {
          n->temp = otmp; // Function output is a duplicate, use the node encountered first
        }

        // Not in the algorithm
        place_in_alg.push_back(-1);
      }
    }

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    vector<casadi_int>& place = place_in_alg; // Reuse memory as it is no longer needed
    place.resize(nodes.size());

    // Stack with unused elements in the work vector, sorted by sparsity pattern
    SPARSITY_MAP<casadi_int, stack<casadi_int> > unused_all;

    // Work vector size
    casadi_int worksize = 0;

    // Find a place in the work vector for the operation
    for (auto&& e : algorithm_) {

      // There are two tasks, allocate memory of the result and free the
      // memory off the arguments, order depends on whether inplace is possible
      casadi_int first_to_free = 0;
      casadi_int last_to_free = e.data->n_inplace();
      for (casadi_int task=0; task<2; ++task) {

        // Dereference or free the memory of the arguments
        for (casadi_int c=last_to_free-1; c>=first_to_free; --c) { // reverse order so that the
                                                          // first argument will end up
                                                          // at the top of the stack

          // Index of the argument
          casadi_int& ch_ind = e.arg[c];
          if (ch_ind>=0) {

            // Decrease reference count and add to the stack of
            // unused variables if the count hits zero
            casadi_int remaining = --refcount[ch_ind];

            // Free variable for reuse
            if (live_variables && remaining==0) {

              // Get a pointer to the sparsity pattern of the argument that can be freed
              casadi_int nnz = nodes[ch_ind]->sparsity().nnz();

              // Add to the stack of unused work vector elements for the current sparsity
              unused_all[nnz].push(place[ch_ind]);
            }

            // Point to the place in the work vector instead of to the place in the list of nodes
            ch_ind = place[ch_ind];
          }
        }

        // Nothing more to allocate
        if (task==1) break;

        // Free the rest in the next iteration
        first_to_free = last_to_free;
        last_to_free = e.arg.size();

        // Allocate/reuse memory for the results of the operation
        for (casadi_int c=0; c<e.res.size(); ++c) {
          if (e.res[c]>=0) {

            // Are reuse of variables (live variables) enabled?
            if (live_variables) {
              // Get a pointer to the sparsity pattern node
              casadi_int nnz = e.data->sparsity(c).nnz();

              // Get a reference to the stack for the current sparsity
              stack<casadi_int>& unused = unused_all[nnz];

              // Try to reuse a variable from the stack if possible (last in, first out)
              if (!unused.empty()) {
                e.res[c] = place[e.res[c]] = unused.top();
                unused.pop();
                continue; // Success, no new element needed in the work vector
              }
            }

            // Allocate a new element in the work vector
            e.res[c] = place[e.res[c]] = worksize++;
          }
        }
      }
    }

    if (verbose_) {
      if (live_variables) {
        casadi_message("Using live variables: work array is " + str(worksize)
                       + " instead of " + str(nodes.size()));
      } else {
        casadi_message("Live variables disabled.");
      }
    }

    // Allocate work vectors (numeric)
    workloc_.resize(worksize+1);
    fill(workloc_.begin(), workloc_.end(), -1);
    size_t wind=0, sz_w=0;
    for (auto&& e : algorithm_) {
      if (e.op!=OP_OUTPUT) {
        for (casadi_int c=0; c<e.res.size(); ++c) {
          if (e.res[c]>=0) {
            alloc_arg(e.data->sz_arg());
            alloc_res(e.data->sz_res());
            alloc_iw(e.data->sz_iw());
            sz_w = max(sz_w, e.data->sz_w());
            if (workloc_[e.res[c]] < 0) {
              workloc_[e.res[c]] = wind;
              wind += e.data->sparsity(c).nnz();
            }
          }
        }
      }
    }
    workloc_.back()=wind;
    for (casadi_int i=0; i<workloc_.size(); ++i) {
      if (workloc_[i]<0) workloc_[i] = i==0 ? 0 : workloc_[i-1];
      workloc_[i] += sz_w;
    }
    sz_w += wind;
    alloc_w(sz_w);

    // Reset the temporary variables
    for (casadi_int i=0; i<nodes.size(); ++i) {
      if (nodes[i]) {
        nodes[i]->temp = 0;
      }
    }

    // Now mark each input's place in the algorithm
    for (auto it=symb_loc.begin(); it!=symb_loc.end(); ++it) {
      it->second->temp = it->first+1;
    }

    // Add input instructions, loop over inputs
    for (casadi_int ind=0; ind<in_.size(); ++ind) {
      // Loop over symbolic primitives of each input
      vector<MX> prim = in_[ind].primitives();
      casadi_int nz_offset=0;
      for (casadi_int p=0; p<prim.size(); ++p) {
        casadi_int i = prim[p].get_temp()-1;
        if (i>=0) {
          // Mark read
          prim[p].set_temp(0);

          // Replace parameter with input instruction
          algorithm_[i].data.own(new Input(prim[p].sparsity(), ind, p, nz_offset));
          algorithm_[i].op = OP_INPUT;
        }
        nz_offset += prim[p]->nnz();
      }
    }

    // Locate free variables
    free_vars_.clear();
    for (auto it=symb_loc.begin(); it!=symb_loc.end(); ++it) {
      casadi_int i = it->second->temp-1;
      if (i>=0) {
        // Save to list of free parameters
        free_vars_.push_back(MX::create(it->second));

        // Remove marker
        it->second->temp=0;
      }
    }

    // Does any embedded function have reference counting for codegen?
    for (auto&& a : algorithm_) {
      if (a.data->has_refcount()) {
        has_refcount_ = true;
        break;
      }
    }
  }

  int MXFunction::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::eval");
    // Work vector and temporaries to hold pointers to operation input and outputs
    const double** arg1 = arg+n_in_;
    double** res1 = res+n_out_;

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      std::stringstream ss;
      disp(ss, false);
      casadi_error("Cannot evaluate \"" + ss.str() + "\" since variables "
                   + str(free_vars_) + " are free.");
    }

    // Evaluate all of the nodes of the algorithm:
    // should only evaluate nodes that have not yet been calculated!
    for (auto&& e : algorithm_) {
      if (e.op==OP_INPUT) {
        // Pass an input
        double *w1 = w+workloc_[e.res.front()];
        casadi_int nnz=e.data.nnz();
        casadi_int i=e.data->ind();
        casadi_int nz_offset=e.data->offset();
        if (arg[i]==nullptr) {
          fill(w1, w1+nnz, 0);
        } else {
          copy(arg[i]+nz_offset, arg[i]+nz_offset+nnz, w1);
        }
      } else if (e.op==OP_OUTPUT) {
        // Get an output
        double *w1 = w+workloc_[e.arg.front()];
        casadi_int nnz=e.data->dep().nnz();
        casadi_int i=e.data->ind();
        casadi_int nz_offset=e.data->offset();
        if (res[i]) copy(w1, w1+nnz, res[i]+nz_offset);
      } else {
        // Point pointers to the data corresponding to the element
        for (casadi_int i=0; i<e.arg.size(); ++i)
          arg1[i] = e.arg[i]>=0 ? w+workloc_[e.arg[i]] : nullptr;
        for (casadi_int i=0; i<e.res.size(); ++i)
          res1[i] = e.res[i]>=0 ? w+workloc_[e.res[i]] : nullptr;

        // Evaluate
        if (e.data->eval(arg1, res1, iw, w)) return 1;
      }
    }
    return 0;
  }

  string MXFunction::print(const AlgEl& el) const {
    stringstream s;
    if (el.op==OP_OUTPUT) {
      s << "output[" << el.data->ind() << "][" << el.data->segment() << "]"
        << " = @" << el.arg.at(0);
    } else if (el.op==OP_SETNONZEROS || el.op==OP_ADDNONZEROS) {
      if (el.res.front()!=el.arg.at(0)) {
        s << "@" << el.res.front() << " = @" << el.arg.at(0) << "; ";
      }
      vector<string> arg(2);
      arg[0] = "@" + str(el.res.front());
      arg[1] = "@" + str(el.arg.at(1));
      s << el.data->disp(arg);
    } else {
      if (el.res.size()==1) {
        s << "@" << el.res.front() << " = ";
      } else {
        s << "{";
        for (casadi_int i=0; i<el.res.size(); ++i) {
          if (i!=0) s << ", ";
          if (el.res[i]>=0) {
            s << "@" << el.res[i];
          } else {
            s << "NULL";
          }
        }
        s << "} = ";
      }
      vector<string> arg;
      if (el.op!=OP_INPUT) {
        arg.resize(el.arg.size());
        for (casadi_int i=0; i<el.arg.size(); ++i) {
          if (el.arg[i]>=0) {
            arg[i] = "@" + str(el.arg[i]);
          } else {
            arg[i] = "NULL";
          }
        }
      }
      s << el.data->disp(arg);
    }
    return s.str();
  }

  void MXFunction::disp_more(ostream &stream) const {
    stream << "Algorithm:";
    for (auto&& e : algorithm_) {
      InterruptHandler::check();
      stream << endl << print(e);
    }
  }

  int MXFunction::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    // Temporaries to hold pointers to operation input and outputs
    const bvec_t** arg1=arg+n_in_;
    bvec_t** res1=res+n_out_;

    // Propagate sparsity forward
    for (auto&& e : algorithm_) {
      if (e.op==OP_INPUT) {
        // Pass input seeds
        casadi_int nnz=e.data.nnz();
        casadi_int i=e.data->ind();
        casadi_int nz_offset=e.data->offset();
        const bvec_t* argi = arg[i];
        bvec_t* w1 = w + workloc_[e.res.front()];
        if (argi!=nullptr) {
          copy(argi+nz_offset, argi+nz_offset+nnz, w1);
        } else {
          fill_n(w1, nnz, 0);
        }
      } else if (e.op==OP_OUTPUT) {
        // Get the output sensitivities
        casadi_int nnz=e.data.dep().nnz();
        casadi_int i=e.data->ind();
        casadi_int nz_offset=e.data->offset();
        bvec_t* resi = res[i];
        bvec_t* w1 = w + workloc_[e.arg.front()];
        if (resi!=nullptr) copy(w1, w1+nnz, resi+nz_offset);
      } else {
        // Point pointers to the data corresponding to the element
        for (casadi_int i=0; i<e.arg.size(); ++i)
          arg1[i] = e.arg[i]>=0 ? w+workloc_[e.arg[i]] : nullptr;
        for (casadi_int i=0; i<e.res.size(); ++i)
          res1[i] = e.res[i]>=0 ? w+workloc_[e.res[i]] : nullptr;

        // Propagate sparsity forwards
        if (e.data->sp_forward(arg1, res1, iw, w)) return 1;
      }
    }
    return 0;
  }

  int MXFunction::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    // Temporaries to hold pointers to operation input and outputs
    bvec_t** arg1=arg+n_in_;
    bvec_t** res1=res+n_out_;

    fill_n(w, sz_w(), 0);

    // Propagate sparsity backwards
    for (auto it=algorithm_.rbegin(); it!=algorithm_.rend(); it++) {
      if (it->op==OP_INPUT) {
        // Get the input sensitivities and clear it from the work vector
        casadi_int nnz=it->data.nnz();
        casadi_int i=it->data->ind();
        casadi_int nz_offset=it->data->offset();
        bvec_t* argi = arg[i];
        bvec_t* w1 = w + workloc_[it->res.front()];
        if (argi!=nullptr) for (casadi_int k=0; k<nnz; ++k) argi[nz_offset+k] |= w1[k];
        fill_n(w1, nnz, 0);
      } else if (it->op==OP_OUTPUT) {
        // Pass output seeds
        casadi_int nnz=it->data.dep().nnz();
        casadi_int i=it->data->ind();
        casadi_int nz_offset=it->data->offset();
        bvec_t* resi = res[i] ? res[i] + nz_offset : nullptr;
        bvec_t* w1 = w + workloc_[it->arg.front()];
        if (resi!=nullptr) {
          for (casadi_int k=0; k<nnz; ++k) w1[k] |= resi[k];
          fill_n(resi, nnz, 0);

        }
      } else {
        // Point pointers to the data corresponding to the element
        for (casadi_int i=0; i<it->arg.size(); ++i)
          arg1[i] = it->arg[i]>=0 ? w+workloc_[it->arg[i]] : nullptr;
        for (casadi_int i=0; i<it->res.size(); ++i)
          res1[i] = it->res[i]>=0 ? w+workloc_[it->res[i]] : nullptr;

        // Propagate sparsity backwards
        if (it->data->sp_reverse(arg1, res1, iw, w)) return 1;
      }
    }
    return 0;
  }

  std::vector<MX> MXFunction::symbolic_output(const std::vector<MX>& arg) const {
    // Check if input is given
    const casadi_int checking_depth = 2;
    bool input_given = true;
    for (casadi_int i=0; i<arg.size() && input_given; ++i) {
      if (!is_equal(arg[i], in_[i], checking_depth)) {
        input_given = false;
      }
    }

    // Return output if possible, else fall back to base class
    if (input_given) {
      return out_;
    } else {
      return FunctionInternal::symbolic_output(arg);
    }
  }

  void MXFunction::eval_mx(const MXVector& arg, MXVector& res,
                           bool always_inline, bool never_inline) const {
    if (verbose_) casadi_message(name_ + "::eval_mx");
    try {
      // Resize the number of outputs
      casadi_assert(arg.size()==n_in_, "Wrong number of input arguments");
      res.resize(out_.size());

      // Trivial inline by default if output known
      if (!never_inline && isInput(arg)) {
        copy(out_.begin(), out_.end(), res.begin());
        return;
      }

      // non-inlining call is implemented in the base-class
      if (!should_inline(always_inline, never_inline)) {
        return FunctionInternal::eval_mx(arg, res, false, true);
      }

      // Symbolic work, non-differentiated
      vector<MX> swork(workloc_.size()-1);
      if (verbose_) casadi_message("Allocated work vector");

      // Split up inputs analogous to symbolic primitives
      vector<vector<MX> > arg_split(in_.size());
      for (casadi_int i=0; i<in_.size(); ++i) arg_split[i] = in_[i].split_primitives(arg[i]);

      // Allocate storage for split outputs
      vector<vector<MX> > res_split(out_.size());
      for (casadi_int i=0; i<out_.size(); ++i) res_split[i].resize(out_[i].n_primitives());

      vector<MX> arg1, res1;

      // Loop over computational nodes in forward order
      casadi_int alg_counter = 0;
      for (auto it=algorithm_.begin(); it!=algorithm_.end(); ++it, ++alg_counter) {
        if (it->op == OP_INPUT) {
          swork[it->res.front()] = project(arg_split.at(it->data->ind()).at(it->data->segment()),
                                           it->data.sparsity(), true);
        } else if (it->op==OP_OUTPUT) {
          // Collect the results
          res_split.at(it->data->ind()).at(it->data->segment()) = swork[it->arg.front()];
        } else if (it->op==OP_PARAMETER) {
          // Fetch parameter
          swork[it->res.front()] = it->data;
        } else {
          // Arguments of the operation
          arg1.resize(it->arg.size());
          for (casadi_int i=0; i<arg1.size(); ++i) {
            casadi_int el = it->arg[i]; // index of the argument
            arg1[i] = el<0 ? MX(it->data->dep(i).size()) : swork[el];
          }

          // Perform the operation
          res1.resize(it->res.size());
          it->data->eval_mx(arg1, res1);

          // Get the result
          for (casadi_int i=0; i<res1.size(); ++i) {
            casadi_int el = it->res[i]; // index of the output
            if (el>=0) swork[el] = res1[i];
          }
        }
      }

      // Join split outputs
      for (casadi_int i=0; i<res.size(); ++i) res[i] = out_[i].join_primitives(res_split[i]);
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("eval_mx", e.what());
    }
  }

  void MXFunction::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                std::vector<std::vector<MX> >& fsens) const {
    if (verbose_) casadi_message(name_ + "::ad_forward(" + str(fseed.size())+ ")");
    try {
      // Allocate results
      casadi_int nfwd = fseed.size();
      fsens.resize(nfwd);
      for (casadi_int d=0; d<nfwd; ++d) {
        fsens[d].resize(n_out_);
      }

      // Quick return if no directions
      if (nfwd==0) return;

      // Check if seeds need to have dimensions corrected
      casadi_int npar = 1;
      for (auto&& r : fseed) {
        if (!matching_arg(r, npar)) {
          casadi_assert_dev(npar==1);
          return ad_forward(replace_fseed(fseed, npar), fsens);
        }
      }

      // Check if there are any zero seeds
      for (auto&& r : fseed) {
        if (purgable(r)) {
          // New argument without all-zero directions
          std::vector<std::vector<MX> > fseed_purged, fsens_purged;
          fseed_purged.reserve(nfwd);
          vector<casadi_int> index_purged;
          for (casadi_int d=0; d<nfwd; ++d) {
            if (purgable(fseed[d])) {
              for (casadi_int i=0; i<fsens[d].size(); ++i) {
                fsens[d][i] = MX(size_out(i));
              }
            } else {
              fseed_purged.push_back(fsens[d]);
              index_purged.push_back(d);
            }
          }

          // Call recursively
          ad_forward(fseed_purged, fsens_purged);

          // Fetch result
          for (casadi_int d=0; d<fseed_purged.size(); ++d) {
            fsens[index_purged[d]] = fsens_purged[d];
          }
          return;
        }
      }

      // Work vector, forward derivatives
      std::vector<std::vector<MX> > dwork(workloc_.size()-1);
      fill(dwork.begin(), dwork.end(), std::vector<MX>(nfwd));
      if (verbose_) casadi_message("Allocated derivative work vector (forward mode)");

      // Split up fseed analogous to symbolic primitives
      vector<vector<vector<MX> > > fseed_split(nfwd);
      for (casadi_int d=0; d<nfwd; ++d) {
        fseed_split[d].resize(fseed[d].size());
        for (casadi_int i=0; i<fseed[d].size(); ++i) {
          fseed_split[d][i] = in_[i].split_primitives(fseed[d][i]);
        }
      }

      // Allocate splited forward sensitivities
      vector<vector<vector<MX> > > fsens_split(nfwd);
      for (casadi_int d=0; d<nfwd; ++d) {
        fsens_split[d].resize(out_.size());
        for (casadi_int i=0; i<out_.size(); ++i) {
          fsens_split[d][i].resize(out_[i].n_primitives());
        }
      }

      // Pointers to the arguments of the current operation
      vector<vector<MX> > oseed, osens;
      oseed.reserve(nfwd);
      osens.reserve(nfwd);
      vector<bool> skip(nfwd, false);

      // Loop over computational nodes in forward order
      for (auto&& e : algorithm_) {
        if (e.op == OP_INPUT) {
          // Fetch forward seed
          for (casadi_int d=0; d<nfwd; ++d) {
            dwork[e.res.front()][d] =
              project(fseed_split[d].at(e.data->ind()).at(e.data->segment()),
                                        e.data.sparsity(), true);
          }
        } else if (e.op==OP_OUTPUT) {
          // Collect forward sensitivity
          for (casadi_int d=0; d<nfwd; ++d) {
            fsens_split[d][e.data->ind()][e.data->segment()] = dwork[e.arg.front()][d];
          }
        } else if (e.op==OP_PARAMETER) {
          // Fetch parameter
          for (casadi_int d=0; d<nfwd; ++d) {
            dwork[e.res.front()][d] = MX();
          }
        } else {
          // Get seeds, ignoring all-zero directions
          oseed.clear();
          for (casadi_int d=0; d<nfwd; ++d) {
            // Collect seeds, skipping directions with only zeros
            vector<MX> seed(e.arg.size());
            skip[d] = true; // All seeds are zero?
            for (casadi_int i=0; i<e.arg.size(); ++i) {
              casadi_int el = e.arg[i];
              if (el<0 || dwork[el][d].is_empty(true)) {
                seed[i] = MX(e.data->dep(i).size());
              } else {
                seed[i] = dwork[el][d];
              }
              if (skip[d] && !seed[i].is_zero()) skip[d] = false;
            }
            if (!skip[d]) oseed.push_back(seed);
          }

          // Perform the operation
          osens.resize(oseed.size());
          if (!osens.empty()) {
            fill(osens.begin(), osens.end(), vector<MX>(e.res.size()));
            e.data.ad_forward(oseed, osens);
          }

          // Store sensitivities
          casadi_int d1=0;
          for (casadi_int d=0; d<nfwd; ++d) {
            for (casadi_int i=0; i<e.res.size(); ++i) {
              casadi_int el = e.res[i];
              if (el>=0) {
                dwork[el][d] = skip[d] ? MX(e.data->sparsity(i).size()) : osens[d1][i];
              }
            }
            if (!skip[d]) d1++;
          }
        }
      }

      // Get forward sensitivities
      for (casadi_int d=0; d<nfwd; ++d) {
        for (casadi_int i=0; i<out_.size(); ++i) {
          fsens[d][i] = out_[i].join_primitives(fsens_split[d][i]);
        }
      }
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("ad_forward", e.what());
    }
  }

  void MXFunction::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                std::vector<std::vector<MX> >& asens) const {
    if (verbose_) casadi_message(name_ + "::ad_reverse(" + str(aseed.size())+ ")");
    try {

      // Allocate results
      casadi_int nadj = aseed.size();
      asens.resize(nadj);
      for (casadi_int d=0; d<nadj; ++d) {
        asens[d].resize(n_in_);
      }

      // Quick return if no directions
      if (nadj==0) return;

      // Check if seeds need to have dimensions corrected
      casadi_int npar = 1;
      for (auto&& r : aseed) {
        if (!matching_res(r, npar)) {
          casadi_assert_dev(npar==1);
          return ad_reverse(replace_aseed(aseed, npar), asens);
        }
      }

      // Check if there are any zero seeds
      for (auto&& r : aseed) {
        // If any direction can be skipped
        if (purgable(r)) {
          // New argument without all-zero directions
          std::vector<std::vector<MX> > aseed_purged, asens_purged;
          aseed_purged.reserve(nadj);
          vector<casadi_int> index_purged;
          for (casadi_int d=0; d<nadj; ++d) {
            if (purgable(aseed[d])) {
              for (casadi_int i=0; i<asens[d].size(); ++i) {
                asens[d][i] = MX(size_in(i));
              }
            } else {
              aseed_purged.push_back(asens[d]);
              index_purged.push_back(d);
            }
          }

          // Call recursively
          ad_reverse(aseed_purged, asens_purged);

          // Fetch result
          for (casadi_int d=0; d<aseed_purged.size(); ++d) {
            asens[index_purged[d]] = asens_purged[d];
          }
          return;
        }
      }

      // Split up aseed analogous to symbolic primitives
      vector<vector<vector<MX> > > aseed_split(nadj);
      for (casadi_int d=0; d<nadj; ++d) {
        aseed_split[d].resize(out_.size());
        for (casadi_int i=0; i<out_.size(); ++i) {
          aseed_split[d][i] = out_[i].split_primitives(aseed[d][i]);
        }
      }

      // Allocate splited adjoint sensitivities
      vector<vector<vector<MX> > > asens_split(nadj);
      for (casadi_int d=0; d<nadj; ++d) {
        asens_split[d].resize(in_.size());
        for (casadi_int i=0; i<in_.size(); ++i) {
          asens_split[d][i].resize(in_[i].n_primitives());
        }
      }

      // Pointers to the arguments of the current operation
      vector<vector<MX> > oseed, osens;
      oseed.reserve(nadj);
      osens.reserve(nadj);
      vector<bool> skip(nadj, false);

      // Work vector, adjoint derivatives
      std::vector<std::vector<MX> > dwork(workloc_.size()-1);
      fill(dwork.begin(), dwork.end(), std::vector<MX>(nadj));

      // Loop over computational nodes in reverse order
      for (auto it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
        if (it->op == OP_INPUT) {
          // Get the adjoint sensitivities
          for (casadi_int d=0; d<nadj; ++d) {
            asens_split[d].at(it->data->ind()).at(it->data->segment()) = dwork[it->res.front()][d];
            dwork[it->res.front()][d] = MX();
          }
        } else if (it->op==OP_OUTPUT) {
          // Pass the adjoint seeds
          for (casadi_int d=0; d<nadj; ++d) {
            MX a = project(aseed_split[d].at(it->data->ind()).at(it->data->segment()),
                           it->data.dep().sparsity(), true);
            if (dwork[it->arg.front()][d].is_empty(true)) {
              dwork[it->arg.front()][d] = a;
            } else {
              dwork[it->arg.front()][d] += a;
            }
          }
        } else if (it->op==OP_PARAMETER) {
          // Clear adjoint seeds
          for (casadi_int d=0; d<nadj; ++d) {
            dwork[it->res.front()][d] = MX();
          }
        } else {
          // Collect and reset seeds
          oseed.clear();
          for (casadi_int d=0; d<nadj; ++d) {
            // Can the direction be skipped completely?
            skip[d] = true;

            // Seeds for direction d
            vector<MX> seed(it->res.size());
            for (casadi_int i=0; i<it->res.size(); ++i) {
              // Get and clear seed
              casadi_int el = it->res[i];
              if (el>=0) {
                seed[i] = dwork[el][d];
                dwork[el][d] = MX();
              } else {
                seed[i] = MX();
              }

              // If first time encountered, reset to zero of right dimension
              if (seed[i].is_empty(true)) seed[i] = MX(it->data->sparsity(i).size());

              // If nonzero seeds, keep direction
              if (skip[d] && !seed[i].is_zero()) skip[d] = false;
            }
            // Add to list of derivatives
            if (!skip[d]) oseed.push_back(seed);
          }

          // Get values of sensitivities before addition
          osens.resize(oseed.size());
          casadi_int d1=0;
          for (casadi_int d=0; d<nadj; ++d) {
            if (skip[d]) continue;
            osens[d1].resize(it->arg.size());
            for (casadi_int i=0; i<it->arg.size(); ++i) {
              // Pass seed and reset to avoid counting twice
              casadi_int el = it->arg[i];
              if (el>=0) {
                osens[d1][i] = dwork[el][d];
                dwork[el][d] = MX();
              } else {
                osens[d1][i] = MX();
              }

              // If first time encountered, reset to zero of right dimension
              if (osens[d1][i].is_empty(true)) osens[d1][i] = MX(it->data->dep(i).size());
            }
            d1++;
          }

          // Perform the operation
          if (!osens.empty()) {
            it->data.ad_reverse(oseed, osens);
          }

          // Store sensitivities
          d1=0;
          for (casadi_int d=0; d<nadj; ++d) {
            if (skip[d]) continue;
            for (casadi_int i=0; i<it->arg.size(); ++i) {
              casadi_int el = it->arg[i];
              if (el>=0) {
                if (dwork[el][d].is_empty(true)) {
                  dwork[el][d] = osens[d1][i];
                } else {
                  dwork[el][d] += osens[d1][i];
                }
              }
            }
            d1++;
          }
        }
      }

      // Get adjoint sensitivities
      for (casadi_int d=0; d<nadj; ++d) {
        for (casadi_int i=0; i<in_.size(); ++i) {
          asens[d][i] = in_[i].join_primitives(asens_split[d][i]);
        }
      }
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("ad_reverse", e.what());
    }
  }

  int MXFunction::eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const {
    // Work vector and temporaries to hold pointers to operation input and outputs
    vector<const SXElem*> argp(sz_arg());
    vector<SXElem*> resp(sz_res());

    // Evaluate all of the nodes of the algorithm:
    // should only evaluate nodes that have not yet been calculated!
    for (auto&& a : algorithm_) {
      if (a.op==OP_INPUT) {
        // Pass an input
        SXElem *w1 = w+workloc_[a.res.front()];
        casadi_int nnz=a.data.nnz();
        casadi_int i=a.data->ind();
        casadi_int nz_offset=a.data->offset();
        if (arg[i]==nullptr) {
          std::fill(w1, w1+nnz, 0);
        } else {
          std::copy(arg[i]+nz_offset, arg[i]+nz_offset+nnz, w1);
        }
      } else if (a.op==OP_OUTPUT) {
        // Get the outputs
        SXElem *w1 = w+workloc_[a.arg.front()];
        casadi_int nnz=a.data.dep().nnz();
        casadi_int i=a.data->ind();
        casadi_int nz_offset=a.data->offset();
        if (res[i]) std::copy(w1, w1+nnz, res[i]+nz_offset);
      } else if (a.op==OP_PARAMETER) {
        continue; // FIXME
      } else {
        // Point pointers to the data corresponding to the element
        for (casadi_int i=0; i<a.arg.size(); ++i)
          argp[i] = a.arg[i]>=0 ? w+workloc_[a.arg[i]] : nullptr;
        for (casadi_int i=0; i<a.res.size(); ++i)
          resp[i] = a.res[i]>=0 ? w+workloc_[a.res[i]] : nullptr;

        // Evaluate
        if (a.data->eval_sx(get_ptr(argp), get_ptr(resp), iw, w)) return 1;
      }
    }
    return 0;
  }

  void MXFunction::codegen_declarations(CodeGenerator& g) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation is not possible since variables "
                   + str(free_vars_) + " are free.");
    }

    // Generate code for the embedded functions
    for (auto&& a : algorithm_) {
      a.data->add_dependency(g);
    }
  }

  void MXFunction::codegen_incref(CodeGenerator& g) const {
    set<void*> added;
    for (auto&& a : algorithm_) {
      a.data->codegen_incref(g, added);
    }
  }

  void MXFunction::codegen_decref(CodeGenerator& g) const {
    set<void*> added;
    for (auto&& a : algorithm_) {
      a.data->codegen_decref(g, added);
    }
  }

  void MXFunction::codegen_body(CodeGenerator& g) const {
    // Temporary variables and vectors
    g.init_local("arg1", "arg+" + str(n_in_));
    g.init_local("res1", "res+" + str(n_out_));

    // Declare scalar work vector elements as local variables
    bool first = true;
    for (casadi_int i=0; i<workloc_.size()-1; ++i) {
      casadi_int n=workloc_[i+1]-workloc_[i];
      if (n==0) continue;
      if (first) {
        g << "casadi_real ";
        first = false;
      } else {
        g << ", ";
      }
      /* Could use local variables for small work vector elements here, e.g.:
         ...
         } else if (n<10) {
         g << "w" << i << "[" << n << "]";
         } else {
         ...
      */
      if (!g.codegen_scalars && n==1) {
        g << "w" << i;
      } else {
        g << "*w" << i << "=w+" << workloc_[i];
      }
    }
    if (!first) g << ";\n";

    // Operation number (for printing)
    casadi_int k=0;

    // Names of operation argument and results
    vector<casadi_int> arg, res;

    // Codegen the algorithm
    for (auto&& e : algorithm_) {
      // Generate comment
      if (g.verbose) {
        g << "/* #" << k++ << ": " << print(e) << " */\n";
      }

      // Get the names of the operation arguments
      arg.resize(e.arg.size());
      for (casadi_int i=0; i<e.arg.size(); ++i) {
        casadi_int j=e.arg.at(i);
        if (j>=0 && workloc_.at(j)!=workloc_.at(j+1)) {
          arg.at(i) = j;
        } else {
          arg.at(i) = -1;
        }
      }

      // Get the names of the operation results
      res.resize(e.res.size());
      for (casadi_int i=0; i<e.res.size(); ++i) {
        casadi_int j=e.res.at(i);
        if (j>=0 && workloc_.at(j)!=workloc_.at(j+1)) {
          res.at(i) = j;
        } else {
          res.at(i) = -1;
        }
      }

      // Generate operation
      e.data->generate(g, arg, res);
    }
  }

  void MXFunction::generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const {
    vector<MX> swork(workloc_.size()-1);

    vector<MX> arg1, res1;

    // Get input primitives
    vector<vector<MX> > in_split(in_.size());
    for (casadi_int i=0; i<in_.size(); ++i) in_split[i] = in_[i].primitives();

    // Definition of intermediate variables
    vector<MX> y;
    vector<MX> g;
    vector<vector<MX> > f_G(out_.size());
    for (casadi_int i=0; i<out_.size(); ++i) f_G[i].resize(out_[i].n_primitives());

    // Initial guess for intermediate variables
    vector<MX> x_init;

    // Temporary stringstream
    stringstream ss;

    for (casadi_int algNo=0; algNo<2; ++algNo) {
      for (auto&& e : algorithm_) {
        switch (e.op) {
        case OP_LIFT:
          {
            MX& arg = swork[e.arg.at(0)];
            MX& arg_init = swork[e.arg.at(1)];
            MX& res = swork[e.res.front()];
            switch (algNo) {
            case 0:
              ss.str(string());
              ss << "y" << y.size();
              y.push_back(MX::sym(ss.str(), arg.sparsity()));
              g.push_back(arg);
              res = y.back();
              break;
            case 1:
              x_init.push_back(arg_init);
              res = arg_init;
              break;
            }
            break;
          }
        case OP_INPUT:
          swork[e.res.front()] = in_split.at(e.data->ind()).at(e.data->segment());
          break;
        case OP_PARAMETER:
          swork[e.res.front()] = e.data;
          break;
        case OP_OUTPUT:
          if (algNo==0) {
            f_G.at(e.data->ind()).at(e.data->segment()) = swork[e.arg.front()];
          }
          break;
        default:
          {
            // Arguments of the operation
            arg1.resize(e.arg.size());
            for (casadi_int i=0; i<arg1.size(); ++i) {
              casadi_int el = e.arg[i]; // index of the argument
              arg1[i] = el<0 ? MX(e.data->dep(i).size()) : swork[el];
            }

            // Perform the operation
            res1.resize(e.res.size());
            e.data.eval_mx(arg1, res1);

            // Get the result
            for (casadi_int i=0; i<res1.size(); ++i) {
              casadi_int el = e.res[i]; // index of the output
              if (el>=0) swork[el] = res1[i];
            }
          }
        }
      }
    }

    // Definition of intermediate variables
    vector<MX> f_in = in_;
    f_in.insert(f_in.end(), y.begin(), y.end());
    vector<MX> f_out;
    for (casadi_int i=0; i<out_.size(); ++i) f_out.push_back(out_[i].join_primitives(f_G[i]));
    f_out.insert(f_out.end(), g.begin(), g.end());
    vdef_fcn = Function("lifting_variable_definition", f_in, f_out);

    // Initial guess of intermediate variables
    f_in = in_;
    f_out = x_init;
    vinit_fcn = Function("lifting_variable_guess", f_in, f_out);
  }

  const MX MXFunction::mx_in(casadi_int ind) const {
    return in_.at(ind);
  }

  const std::vector<MX> MXFunction::mx_in() const {
    return in_;
  }

  bool MXFunction::is_a(const std::string& type, bool recursive) const {
    return type=="MXFunction"
      || (recursive && XFunction<MXFunction,
          MX, MXNode>::is_a(type, recursive));
  }

  void MXFunction::substitute_inplace(std::vector<MX>& vdef, std::vector<MX>& ex) const {
    vector<MX> work(workloc_.size()-1);
    vector<MX> oarg, ores;

    for (auto it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        casadi_assert(it->data->segment()==0, "Not implemented");
        work.at(it->res.front()) = vdef.at(it->data->ind());
        break;
      case OP_PARAMETER:
      case OP_CONST:
        work.at(it->res.front()) = it->data;
        break;
      case OP_OUTPUT:
        casadi_assert(it->data->segment()==0, "Not implemented");
        if (it->data->ind()<vdef.size()) {
          vdef.at(it->data->ind()) = work.at(it->arg.front());
        } else {
          ex.at(it->data->ind()-vdef.size()) = work.at(it->arg.front());
        }
        break;
      default:
        {
          // Arguments of the operation
          oarg.resize(it->arg.size());
          for (casadi_int i=0; i<oarg.size(); ++i) {
            casadi_int el = it->arg[i];
            oarg[i] = el<0 ? MX(it->data->dep(i).size()) : work.at(el);
          }

          // Perform the operation
          ores.resize(it->res.size());
          it->data->eval_mx(oarg, ores);

          // Get the result
          for (casadi_int i=0; i<ores.size(); ++i) {
            casadi_int el = it->res[i];
            if (el>=0) work.at(el) = ores[i];
          }
        }
      }
    }
  }

  bool MXFunction::should_inline(bool always_inline, bool never_inline) const {
    // If inlining has been specified
    casadi_assert(!(always_inline && never_inline),
      "Inconsistent options for " + definition());
    casadi_assert(!(never_inline && has_free()),
      "Must inline " + definition());
    if (always_inline) return true;
    if (never_inline) return false;
    // Functions with free variables must be inlined
    if (has_free()) return true;
    // No inlining by default
    return false;
  }

  void MXFunction::export_code_body(const std::string& lang,
      std::ostream &ss, const Dict& options) const {

    // Default values for options
    casadi_int indent_level = 0;

    // Read options
    for (auto&& op : options) {
      if (op.first=="indent_level") {
        indent_level = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Construct indent string
    std::string indent = "";
    for (casadi_int i=0;i<indent_level;++i) {
      indent += "  ";
    }

    Function f = shared_from_this<Function>();

    // Loop over algorithm
    for (casadi_int k=0;k<f.n_instructions();++k) {
      // Get operation
      casadi_int op = static_cast<casadi_int>(f.instruction_id(k));
      // Get MX node
      MX x = f.instruction_MX(k);
      // Get input positions into workvector
      std::vector<casadi_int> o = f.instruction_output(k);
      // Get output positions into workvector
      std::vector<casadi_int> i = f.instruction_input(k);

      switch (op) {
        case OP_INPUT:
          ss << indent << "w" << o[0] << " = varargin{" << i[0]+1 << "};" << std::endl;
          break;
        case OP_OUTPUT:
          {
            Dict info = x.info();
            casadi_int segment = info["segment"];
            x.dep(0).sparsity().export_code("matlab", ss,
              {{"name", "sp_in"}, {"indent_level", indent_level}, {"as_matrix", true}});
            ss << indent << "argout_" << o[0] << "{" << (1+segment) << "} = ";
            ss << "w" << i[0] << "(sp_in==1);" << std::endl;
          }
          break;
        case OP_CONST:
          {
            DM v = static_cast<DM>(x);
            Dict opts;
            opts["name"] = "m";
            opts["indent_level"] = indent_level;
            v.export_code("matlab", ss, opts);
            ss << indent << "w" << o[0] << " = m;" << std::endl;
          }
          break;
        case OP_SQ:
          ss << indent << "w" << o[0] << " = " << "w" << i[0] << ".^2;" << std::endl;
          break;
        case OP_MTIMES:
          ss << indent << "w" << o[0] << " = ";
          ss << "w" << i[1] << "*w" << i[2] << "+w" << i[0] << ";" << std::endl;
          break;
        case OP_MUL:
          {
            std::string prefix = (x.dep(0).is_scalar() || x.dep(1).is_scalar()) ? "" : ".";
            ss << indent << "w" << o[0] << " = " << "w" << i[0] << prefix << "*w" << i[1] << ";";
            ss << std::endl;
          }
          break;
        case OP_TWICE:
          ss << indent << "w" << o[0] << " = 2*w" << i[0] << ";" << std::endl;
          break;
        case OP_INV:
          ss << indent << "w" << o[0] << " = 1./w" << i[0] << ";" << std::endl;
          break;
        case OP_DOT:
          ss << indent << "w" << o[0] << " = dot(w" << i[0] << ",w" << i[1]<< ");" << std::endl;
          break;
        case OP_BILIN:
          ss << indent << "w" << o[0] << " = w" << i[1] << ".'*w" << i[0]<< "*w" << i[2] << ";";
          ss << std::endl;
          break;
        case OP_RANK1:
          ss << indent << "w" << o[0] << " = w" << i[0] << "+";
          ss << "w" << i[1] << "*w" << i[2] << "*w" << i[3] << ".';";
          ss << std::endl;
          break;
        case OP_FABS:
          ss << indent << "w" << o[0] << " = abs(w" << i[0] << ");" << std::endl;
          break;
        case OP_DETERMINANT:
          ss << indent << "w" << o[0] << " = det(w" << i[0] << ");" << std::endl;
          break;
        case OP_INVERSE:
          ss << indent << "w" << o[0] << " = inv(w" << i[0] << ");";
          ss << "w" << o[0] << "(w" << o[0] << "==0) = 1e-200;" << std::endl;
          break;
        case OP_SOLVE:
          {
            bool tr = x.info()["tr"];
            if (tr) {
              ss << indent << "w" << o[0] << " = ((w" << i[1] << ".')\\w" << i[0] << ").';";
              ss << std::endl;
            } else {
              ss << indent << "w" << o[0] << " = w" << i[1] << "\\w" << i[0] << ";" << std::endl;
            }
            ss << "w" << o[0] << "(w" << o[0] << "==0) = 1e-200;" << std::endl;
          }
          break;
        case OP_DIV:
          {
            std::string prefix = (x.dep(0).is_scalar() || x.dep(1).is_scalar()) ? "" : ".";
            ss << indent << "w" << o[0] << " = " << "w" << i[0] << prefix << "/w" << i[1] << ";";
            ss << std::endl;
          }
          break;
        case OP_POW:
        case OP_CONSTPOW:
          ss << indent << "w" << o[0] << " = " << "w" << i[0] << ".^w" << i[1] << ";" << std::endl;
          break;
        case OP_TRANSPOSE:
          ss << indent << "w" << o[0] << " = " << "w" << i[0] << ".';" << std::endl;
          break;
        case OP_HORZCAT:
        case OP_VERTCAT:
          {
            ss << indent << "w" << o[0] << " = [";
            for (casadi_int e : i) {
              ss << "w" << e << (op==OP_HORZCAT ? " " : ";");
            }
            ss << "];" << std::endl;
          }
          break;
        case OP_DIAGCAT:
          {
            for (casadi_int k=0;k<i.size();++k) {
              x.dep(k).sparsity().export_code("matlab", ss,
                {{"name", "sp_in" + str(k)}, {"indent_level", indent_level}, {"as_matrix", true}});
            }
            ss << indent << "w" << o[0] << " = [";
            for (casadi_int k=0;k<i.size();++k) {
              ss << "w" << i[k] << "(sp_in" << k << "==1);";
            }
            ss << "];" << std::endl;
            Dict opts;
            opts["name"] = "sp";
            opts["indent_level"] = indent_level;
            opts["as_matrix"] = false;
            x.sparsity().export_code("matlab", ss, opts);
            ss << indent << "w" << o[0] << " = ";
            ss << "sparse(sp_i, sp_j, w" << o[0] << ", sp_m, sp_n);" << std::endl;
          }
          break;
        case OP_HORZSPLIT:
        case OP_VERTSPLIT:
          {
            Dict info = x.info();
            std::vector<casadi_int> offset = info["offset"];
            casadi::Function output = info["output"];
            std::vector<Sparsity> sp;
            for (casadi_int i=0;i<output.n_out();i++)
              sp.push_back(output.sparsity_out(i));
            for (casadi_int k=0;k<o.size();++k) {
              if (o[k]==-1) continue;
              x.dep(0).sparsity().export_code("matlab", ss,
                {{"name", "sp_in"}, {"indent_level", indent_level}, {"as_matrix", true}});
              ss << indent << "tmp = w" << i[0]<< "(sp_in==1);" << std::endl;
              Dict opts;
              opts["name"] = "sp";
              opts["indent_level"] = indent_level;
              opts["as_matrix"] = false;
              sp[k].export_code("matlab", ss, opts);
              ss << indent << "w" << o[k] << " = sparse(sp_i, sp_j, ";
              ss << "tmp(" << offset[k]+1 << ":" << offset[k+1] << "), sp_m, sp_n);" << std::endl;
            }
          }
          break;
        case OP_GETNONZEROS:
        case OP_SETNONZEROS:
          {
            Dict info = x.info();

            std::string nonzeros;
            if (info.find("nz")!=info.end()) {
              nonzeros = "1+" + str(info["nz"]);
            } else if (info.find("slice")!=info.end()) {
              Dict s = info["slice"];
              casadi_int start = s["start"];
              casadi_int step = s["step"];
              casadi_int stop = s["stop"];
              nonzeros = str(start+1) + ":" + str(step) + ":" + str(stop);
              nonzeros = "nonzeros(" + nonzeros + ")";
            } else {
              Dict inner = info["inner"];
              Dict outer = info["outer"];
              casadi_int inner_start = inner["start"];
              casadi_int inner_step  = inner["step"];
              casadi_int inner_stop  = inner["stop"];
              casadi_int outer_start = outer["start"];
              casadi_int outer_step  = outer["step"];
              casadi_int outer_stop  = outer["stop"];
              std::string inner_slice = "(" + str(inner_start) + ":" +
                str(inner_step) + ":" + str(inner_stop-1)+")";
              std::string outer_slice = "(" + str(outer_start+1) + ":" +
                str(outer_step) + ":" + str(outer_stop)+")";
              casadi_int N = range(outer_start, outer_stop, outer_step).size();
              casadi_int M = range(inner_start, inner_stop, inner_step).size();
              nonzeros = "repmat("+ inner_slice  +"', 1, " + str(N) + ")+" +
                         "repmat("+ outer_slice  +", " + str(M) + ", 1)";
              nonzeros = "nonzeros(" + nonzeros + ")";
            }

            Dict opts;
            opts["name"] = "sp";
            opts["indent_level"] = indent_level;
            opts["as_matrix"] = false;
            x.sparsity().export_code("matlab", ss, opts);

            if (op==OP_GETNONZEROS) {
              x.dep(0).sparsity().export_code("matlab", ss,
                {{"name", "sp_in"}, {"indent_level", indent_level}, {"as_matrix", true}});
              //ss << indent << "w" << i[0] << "" << std::endl;
              //ss << indent << "size(w" << i[0] << ")" << std::endl;
              ss << indent << "in_flat = w" << i[0] << "(sp_in==1);" << std::endl;
              //ss << indent << "in_flat" << std::endl;
              //ss << indent << "size(in_flat)" << std::endl;
              ss << indent << "w" << o[0] << " = in_flat(" << nonzeros << ");" << std::endl;
            } else {
              x.dep(0).sparsity().export_code("matlab", ss,
                {{"name", "sp_in0"}, {"indent_level", indent_level}, {"as_matrix", true}});
              x.dep(1).sparsity().export_code("matlab", ss,
                {{"name", "sp_in1"}, {"indent_level", indent_level}, {"as_matrix", true}});
              ss << indent << "in_flat = w" << i[1] << "(sp_in1==1);" << std::endl;
              ss << indent << "w" << o[0] << " = w" << i[0] << "(sp_in0==1);" << std::endl;
              ss << indent << "w" << o[0] << "(" << nonzeros << ")  = ";
              if (info["add"]) ss << "w" << o[0] << "(" << nonzeros << ") + ";
              ss << "in_flat;";
            }
            ss << indent << "w" << o[0] << " = ";
            ss << "sparse(sp_i, sp_j, w" << o[0] << ", sp_m, sp_n);" << std::endl;
          }
          break;
        case OP_PROJECT:
          {
            Dict opts;
            opts["name"] = "sp";
            opts["indent_level"] = indent_level;
            x.sparsity().export_code("matlab", ss, opts);
            ss << indent << "w" << o[0] << " = ";
            ss << "sparse(sp_i, sp_j, w" << i[0] << "(sp==1), sp_m, sp_n);" << std::endl;
          }
          break;
        case OP_NORM1:
          ss << indent << "w" << o[0] << " = norm(w" << i[0] << ", 1);" << std::endl;
          break;
        case OP_NORM2:
          ss << indent << "w" << o[0] << " = norm(w" << i[0] << ", 2);" << std::endl;
          break;
        case OP_NORMF:
          ss << indent << "w" << o[0] << " = norm(w" << i[0] << ", 'fro');" << std::endl;
          break;
        case OP_NORMINF:
          ss << indent << "w" << o[0] << " = norm(w" << i[0] << ", inf);" << std::endl;
          break;
        case OP_MMIN:
          ss << indent << "w" << o[0] << " = min(w" << i[0] << ");" << std::endl;
          break;
        case OP_MMAX:
          ss << indent << "w" << o[0] << " = max(w" << i[0] << ");" << std::endl;
          break;
        case OP_NOT:
          ss << indent << "w" << o[0] << " = ~" << "w" << i[0] << ";" << std::endl;
          break;
        case OP_OR:
          ss << indent << "w" << o[0] << " = w" << i[0] << " | w" << i[1] << ";" << std::endl;
          break;
        case OP_AND:
          ss << indent << "w" << o[0] << " = w" << i[0] << " & w" << i[1] << ";" << std::endl;
          break;
        case OP_NE:
          ss << indent << "w" << o[0] << " = w" << i[0] << " ~= w" << i[1] << ";" << std::endl;
          break;
        case OP_IF_ELSE_ZERO:
          ss << indent << "w" << o[0] << " = ";
          ss << "if_else_zero_gen(w" << i[0] << ", w" << i[1] << ");" << std::endl;
          break;
        case OP_RESHAPE:
          {
            x.dep(0).sparsity().export_code("matlab", ss,
              {{"name", "sp_in"}, {"indent_level", indent_level}, {"as_matrix", true}});
            x.sparsity().export_code("matlab", ss,
              {{"name", "sp_out"}, {"indent_level", indent_level}, {"as_matrix", false}});
            ss << indent << "w" << o[0] << " = sparse(sp_out_i, sp_out_j, ";
            ss << "w" << i[0] << "(sp_in==1), sp_out_m, sp_out_n);" << std::endl;
          }
          break;
        default:
          if (x.is_binary()) {
            ss << indent << "w" << o[0] << " = " << casadi::casadi_math<double>::print(op,
              "w"+std::to_string(i[0]), "w"+std::to_string(i[1])) << ";" << std::endl;
          } else if (x.is_unary()) {
            ss << indent << "w" << o[0] << " = " << casadi::casadi_math<double>::print(op,
              "w"+std::to_string(i[0])) << ";" << std::endl;
          } else {
            ss << "unknown" + x.class_name() << std::endl;
          }
      }
    }
  }

  Dict MXFunction::get_stats(void* mem) const {
    Dict stats = XFunction::get_stats(mem);

    Function dep;
    for (auto&& e : algorithm_) {
      if (e.op==OP_CALL) {
        Function d = e.data.which_function();
        if (d.is_a("conic", true)) {
          if (!dep.is_null()) return {};
          dep = d;
        }
      }
    }
    if (dep.is_null()) return {};
    return dep.stats(1);
  }

} // namespace casadi
