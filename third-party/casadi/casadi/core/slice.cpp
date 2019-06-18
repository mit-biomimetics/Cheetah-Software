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


#include "slice.hpp"
#include "casadi_misc.hpp"

using namespace std;
namespace casadi {

  Slice::Slice() : start(0), stop(std::numeric_limits<casadi_int>::max()), step(1) {
  }

  Slice::Slice(casadi_int i, bool ind1) : start(i-ind1), stop(i-ind1+1), step(1) {
    casadi_assert(!(ind1 && i<=0),
                          "Matlab is 1-based, but requested index " +
                          str(i) +  ". Note that negative slices are"
                          " disabled in the Matlab interface. "
                          "Possibly you may want to use 'end'.");
    if (i==-1) stop = std::numeric_limits<casadi_int>::max();
  }

  Slice::Slice(casadi_int start, casadi_int stop, casadi_int step) :
    start(start), stop(stop), step(step) { }

  Slice::Slice(int start, int stop, int step) : start(start), stop(stop), step(step) {
  }
  Slice::Slice(int start, casadi_int stop, int step) : start(start), stop(stop), step(step) {
  }
  Slice::Slice(casadi_int start, int stop, int step) : start(start), stop(stop), step(step) {
  }

  std::vector<casadi_int> Slice::all(casadi_int len, bool ind1) const {
    casadi_int start = this->start;
    if (start==std::numeric_limits<casadi_int>::min()) {
      start = (step < 0) ? len - 1 : 0;
    } else if (start<0) {
      start+=len;
    }
    casadi_int stop = this->stop;
    if (stop==std::numeric_limits<casadi_int>::max()) {
      stop = (step < 0) ? -1 : len;
    } else if (stop<0) {
      stop+=len;
    }

    casadi_assert(stop<=len,
      "Slice (start=" + str(start) + ", stop=" + str(stop) + ", step=" + str(step)
      + ") out of bounds with supplied length of " + str(len));
    casadi_assert(start>=0,
      "Slice (start=" + str(start) + ", stop=" + str(stop) + ", step=" + str(step)
      + ") out of bounds with start<0.");
    if ((stop>=start && step<0) || (stop<=start && step>0)) return std::vector<casadi_int>();

    return range(start+ind1, stop+ind1, step, len+ind1);
  }

  void Slice::disp(std::ostream& stream, bool more) const {
    bool from_beginning = start == 0;
    bool till_end = stop == std::numeric_limits<casadi_int>::max();
    bool skip_none = step==1;
    if (stop==start+1) {
      stream << start;
    } else {
      if (!from_beginning) stream << start;
      stream << ":";
      if (!till_end) stream << stop;
      if (!skip_none) stream << ":" << step;
    }
  }

  std::vector<casadi_int> Slice::all(const Slice& outer, casadi_int len) const {
    std::vector<casadi_int> ret;
    for (casadi_int i=outer.start; i!=outer.stop; i+=outer.step) {
      for (casadi_int j=i+start; j!=i+stop; j+=step) {
        ret.push_back(j);
      }
    }
    return ret;
  }

  bool Slice::is_scalar(casadi_int len) const {
    casadi_int start = std::min(this->start, len);
    casadi_int stop = std::min(this->stop, len);
    casadi_int nret = (stop-start)/step + ((stop-start)%step!=0);
    return nret==1;
  }

  casadi_int Slice::scalar(casadi_int len) const {
    casadi_assert_dev(is_scalar(len));
    casadi_assert(start >= -len && start < len, "Slice::getScalar: out of bounds");
    return start >= 0 ? start : start+len;
  }

  Slice CASADI_EXPORT to_slice(const std::vector<casadi_int>& v, bool ind1) {
    Slice r;
    casadi_assert(is_slice(v, ind1), "Cannot be represented as a Slice");
    if (v.size()==0) {
      r.start=r.stop=0;
      r.step = 1;
    } else if (v.size()==1) {
      r.start = v.front()-ind1;
      r.stop = r.start + 1;
      r.step = 1;
    } else {
      r.start = v[0]-ind1;
      r.step = v[1]-v[0];
      r.stop = r.start + r.step*v.size();
    }
    return r;
  }

  bool CASADI_EXPORT is_slice(const std::vector<casadi_int>& v, bool ind1) {
    // Always false if negative numbers or non-increasing
    casadi_int last_v = -1;
    for (casadi_int i=0; i<v.size(); ++i) {
      casadi_assert(!(ind1 && v[i]<=0),
        "Matlab is 1-based, but requested index " + str(v[i]) + ". "
        "Note that negative slices are disabled in the Matlab interface. "
        "Possibly you may want to use 'end'.");
      if (v[i]-ind1<=last_v) return false;
      last_v = v[i]-ind1;
    }

    // Always true if less than 2 elements
    if (v.size()<2) return true;

    // If two elements, true if they are different
    if (v.size()==2) return v[0]!=v[1];

    // We can now get the beginning, end and step
    casadi_int start = v[0]-ind1;
    casadi_int step = v[1]-v[0];
    //casadi_int stop = start + step*v.size();

    // Consistency check
    for (casadi_int i=2; i<v.size(); ++i) {
      if (v[i]-ind1!=start+i*step) return false;
    }

    // True if reached this point
    return true;
  }

  bool CASADI_EXPORT is_slice2(const std::vector<casadi_int>& v) {
    // Always true if 1D slice
    if (is_slice(v)) return true;

    // Always false if negative numbers or non-increasing
    casadi_int last_v = -1;
    for (casadi_int i=0; i<v.size(); ++i) {
      if (v[i]<=last_v) return false;
      last_v = v[i];
    }

    // Get the slices
    casadi_int start_outer = 0;
    casadi_int step_outer = -1;
    casadi_int start_inner = v.front();
    casadi_int step_inner = v[1]-v[0];
    casadi_int stop_inner = -1;
    for (casadi_int i=2; i<v.size(); ++i) {
      casadi_int predicted_v = start_inner+i*step_inner;
      if (v[i]!=predicted_v) {
        stop_inner = predicted_v;
        step_outer = v[i] - start_inner;
        break;
      }
    }
    casadi_assert_dev(stop_inner>=0);

    // Get the end of the outer slice
    casadi_int stop_outer = v.back();
    do {
      if (step_outer>0) stop_outer++;
      else             stop_outer--;
    } while (stop_outer % step_outer!=0);

    // Check consistency
    std::vector<casadi_int>::const_iterator it=v.begin();
    for (casadi_int i=start_outer; i!=stop_outer; i+=step_outer) {
      for (casadi_int j=i+start_inner; j!=i+stop_inner; j+=step_inner) {
        // False if we've reached the end
        if (it==v.end()) return false;

        // Check if value matches
        if (*it++ != j) return false;
      }
    }

    // False if there are still elements not accounted for
    if (it!=v.end()) return false;

    // True if reached this point
    return true;
  }

  std::pair<Slice, Slice> CASADI_EXPORT to_slice2(const std::vector<casadi_int>& v) {
    casadi_assert(is_slice2(v), "Cannot be represented as a nested Slice");
    Slice inner, outer;

    // If simple slice
    if (is_slice(v)) {
      inner = to_slice(v);
      outer.start = 0;
      outer.step = outer.stop = inner.stop;
      return make_pair(inner, outer);
    }

    // Get the slices
    outer.start = 0;
    outer.step = -1;
    inner.start = v.front();
    inner.step = v[1]-v[0];
    inner.stop = -1;
    for (casadi_int i=2; i<v.size(); ++i) {
      casadi_int predicted_v = inner.start+i*inner.step;
      if (v[i]!=predicted_v) {
        inner.stop = predicted_v;
        outer.step = v[i] - inner.start;
        break;
      }
    }

    // Get the end of the outer slice
    outer.stop = v.back();
    do {
      if (outer.step>0) outer.stop++;
      else              outer.stop--;
    } while (outer.stop % outer.step!=0);
    return make_pair(inner, outer);
  }

} // namespace casadi
