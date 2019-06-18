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


#include "interpolant_impl.hpp"
#include "casadi_misc.hpp"
#include "mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_interpolant(const string& name) {
    return Interpolant::has_plugin(name);
  }

  void load_interpolant(const string& name) {
    Interpolant::load_plugin(name);
  }

  string doc_interpolant(const string& name) {
    return Interpolant::getPlugin(name).doc;
  }

  Function interpolant(const std::string& name,
                       const std::string& solver,
                       const std::vector<std::vector<double> >& grid,
                       const std::vector<double>& values,
                       const Dict& opts) {

    // Dimension at least 1
    casadi_assert(grid.size()>0, "At least one input required");

    // Consistency check, number of elements
    casadi_uint nel=1;
    for (auto&& g : grid) {
      casadi_assert(g.size()>=2, "Need at least two grid points for every input");
      nel *= g.size();
    }
    casadi_assert(values.size() % nel== 0, "Inconsistent number of elements");
    casadi_assert(values.size()>0, "Values cannot be empty");

    // Grid must be strictly increasing
    for (auto&& g : grid) {
      double last = -inf;
      for (auto&& e : g) {
        casadi_assert(!isinf(e) && e>last,
          "Gridpoints must be finite and strictly increasing");
        last = e;
      }
    }

    // Get offset for each input dimension
    vector<casadi_int> offset;
    offset.reserve(grid.size()+1);
    offset.push_back(0);
    for (auto&& g : grid) offset.push_back(offset.back()+g.size());

    // Stack input grids
    vector<double> stacked;
    stacked.reserve(offset.back());
    for (auto&& g : grid) stacked.insert(stacked.end(), g.begin(), g.end());
    casadi_int m = values.size()/nel;
    return Function::create(Interpolant::getPlugin(solver)
                            .creator(name, stacked, offset, values, m), opts);
  }

  Interpolant::
  Interpolant(const std::string& name,
              const std::vector<double>& grid,
              const std::vector<casadi_int>& offset,
              const std::vector<double>& values,
              casadi_int m)
              : FunctionInternal(name), m_(m), grid_(grid), offset_(offset),  values_(values) {
    // Number of grid points
    ndim_ = offset_.size()-1;
  }

  Interpolant::~Interpolant() {
  }

  Sparsity Interpolant::get_sparsity_in(casadi_int i) {
    casadi_assert_dev(i==0);
    return Sparsity::dense(ndim_);
  }

  Sparsity Interpolant::get_sparsity_out(casadi_int i) {
    casadi_assert_dev(i==0);
    return Sparsity::dense(m_);
  }

  std::string Interpolant::get_name_in(casadi_int i) {
    casadi_assert_dev(i==0);
    return "x";
  }

  std::string Interpolant::get_name_out(casadi_int i) {
    casadi_assert_dev(i==0);
    return "f";
  }

  std::map<std::string, Interpolant::Plugin> Interpolant::solvers_;

  const std::string Interpolant::infix_ = "interpolant";

  Options Interpolant::options_
  = {{&FunctionInternal::options_},
     {{"lookup_mode",
       {OT_STRINGVECTOR,
        "Specifies, for each grid dimenion, the lookup algorithm used to find the correct index. "
        "'linear' uses a for-loop + break; (default when #knots<=100), "
        "'exact' uses floored division (only for uniform grids), "
        "'binary' uses a binary search. (default when #knots>100)."}}
     }
  };

  void Interpolant::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_modes_ = op.second;
      }
    }

    // Needed by casadi_interpn
    alloc_w(ndim_, true);
    alloc_iw(2*ndim_, true);
  }

  std::vector<std::string> Interpolant::lookup_mode_from_enum(
      const std::vector<casadi_int>& modes) {
    std::vector<std::string> ret(modes.size());
    for (casadi_int i=0;i<modes.size();++i) {
      switch (modes[i]) {
        case 0:
          ret[i] = "linear";
          break;
        case 1:
          ret[i] = "exact";
          break;
        case 2:
          ret[i] = "binary";
          break;
        default:
          casadi_error("lookup_mode error.");
      }
    }
    return ret;
  }

  std::vector<casadi_int> Interpolant::interpret_lookup_mode(
      const std::vector<std::string>& modes, const std::vector<double>& knots,
      const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& margin_left, const std::vector<casadi_int>& margin_right) {

    // Default lookup mode linear
    std::vector<casadi_int> ret(offset.size()-1, 0);

    for (casadi_int i=0;i<ret.size();++i) {
      // If more than 100 knots -> default is binary search
      if (offset[i+1]-offset[i]>100) ret[i] = 2;
    }

    if (modes.empty()) return ret;

    casadi_assert_dev(modes.size()==offset.size()-1);
    for (casadi_int i=0;i<offset.size()-1;++i) {
      if (modes[i]=="linear") {
        ret[i] = 0;
      } else if (modes[i]=="exact") {
        ret[i] = 1;
        casadi_int m_left  = margin_left.empty() ? 0 : margin_left[i];
        casadi_int m_right = margin_right.empty() ? 0 : margin_right[i];

        std::vector<double> grid(
            knots.begin()+offset[i]+m_left,
            knots.begin()+offset[i+1]-m_right);
        casadi_assert_dev(is_increasing(grid) && is_equally_spaced(grid));
      } else if (modes[i]=="binary") {
        ret[i] = 2;
      } else {
        casadi_error("Unknown lookup_mode option '" + modes[i] + ". "
                     "Allowed values: linear, binary, exact.");
      }
    }
    return ret;
  }
} // namespace casadi
