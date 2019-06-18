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


#include "mx.hpp"

#include "casadi_misc.hpp"
#ifdef HAVE_MKSTEMPS
#include <unistd.h>
#else // HAVE_MKSTEMPS
#ifdef HAVE_SIMPLE_MKSTEMPS
#ifdef _WIN32
#include <io.h>
#include <share.h>
#else
#include <unistd.h>
#endif
#include <random>
#include <chrono>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#endif // HAVE_SIMPLE_MKSTEMPS
#endif // HAVE_MKSTEMPS
using namespace std;

namespace casadi {

  int to_int(casadi_int rhs) {
    casadi_assert(rhs<=std::numeric_limits<int>::max(), "Integer overflow detected.");
    casadi_assert(rhs>=std::numeric_limits<int>::min(), "Integer overflow detected.");
    return rhs;
  }

  std::vector<int> to_int(const std::vector<casadi_int>& rhs) {
    std::vector<int> ret;
    ret.reserve(rhs.size());
    for (casadi_int e : rhs) ret.push_back(to_int(e));
    return ret;
  }

  std::vector< std::vector<int> > to_int(
      const std::vector< std::vector<casadi_int> >& rhs) {
    std::vector< std::vector<int> > ret;
    ret.reserve(rhs.size());
    for (const std::vector<casadi_int>& e : rhs) ret.push_back(to_int(e));
    return ret;
  }

  bool all(const std::vector<bool>& v) {
    for (auto && e : v) {
      if (!e) return false;
    }
    return true;
  }

  bool any(const std::vector<bool>& v) {
    for (auto && e : v) {
      if (e) return true;
    }
    return false;
  }

  std::vector<casadi_int> range(casadi_int start, casadi_int stop,
      casadi_int step, casadi_int len) {
    start = std::min(start, len);
    stop = std::min(stop, len);
    casadi_int nret = (stop-start)/step + ((stop-start)%step!=0);
    std::vector<casadi_int> ret(nret);
    casadi_int ind = start;
    for (std::vector<casadi_int>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = ind;
      ind += step;
    }
    return ret;
  }

  bool is_equally_spaced(const std::vector<double> &v) {
    if (v.size()<=1) return true;

    double margin = (v[v.size()-1]-v[0])*1e-14;

    for (casadi_int i=2;i<v.size();++i) {
      double ref = v[0]+(static_cast<double>(i)*(v[v.size()-1]-v[0]))/
        static_cast<double>(v.size()-1);
      if (abs(ref-v[i])>margin) return false;
    }
    return true;
  }

  std::vector<casadi_int> range(casadi_int stop) {
    return range(0, stop);
  }

  std::vector<casadi_int> complement(const std::vector<casadi_int> &v, casadi_int size) {
    casadi_assert(in_range(v, size),
                          "complement: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<casadi_int> lookup(size, 0);
    std::vector<casadi_int> ret;

    for (casadi_int i=0;i<v.size();i++) {
      lookup[v[i]] = 1;
    }

    for (casadi_int i=0;i<size;i++) {
      if (lookup[i]==0) ret.push_back(i);
    }

    return ret;

  }

  std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v, casadi_int size) {
    casadi_assert(in_range(v, size),
                          "lookupvector: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<casadi_int> lookup(size, -1);

    for (casadi_int i=0;i<v.size();i++) {
      lookup[v[i]] = i;
    }
    return lookup;
  }

  std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v) {
    casadi_assert_dev(!has_negative(v));
    return lookupvector(v, (*std::max_element(v.begin(), v.end()))+1);
  }

  // Better have a bool return flag saying if we need reorer at all
  std::vector<casadi_int> tensor_permute_mapping(const std::vector<casadi_int>& dims,
      const std::vector<casadi_int>& order) {

     // Get problem dimensions
     casadi_int N = casadi::product(dims);
     casadi_int n = dims.size();
     // Quick return if no elements
     if (N==0) return std::vector<casadi_int>();

     // One dimesion => null-permutation
     if (n==1) return range(N);

     // Allocate space for resulting mapping
     std::vector<casadi_int> mapping(N);
     // Quick return if scalar
     if (n==0) return mapping;


     // Compute cumulative product
     std::vector<casadi_int> cumprod(n+1, 1);
     for (casadi_int k=1;k<dims.size();++k) cumprod[k]=cumprod[k-1]*dims[k-1];

     // Elementary stride
     casadi_int stride = cumprod[order[0]];

     // Split problem in inner and outer part
     casadi_int N_inner = dims[order[0]];
     casadi_int N_outer = N/N_inner;

     // Reorder dims, cumprod
     std::vector<casadi_int> new_dims(n-1), new_cumprod(n-1, 1);
     for (casadi_int k=0;k<n-1;++k) {
       new_dims[k] = dims[order[k+1]];
       new_cumprod[k] = cumprod[order[k+1]];
     }

     // Bank of counters
     std::vector<casadi_int> index_counters(n-1);

     // Inex into mapping
     casadi_int m_ind = 0;

     for (casadi_int i=0;i<N_outer;++i) {
       // Compute index
       casadi_int ind = 0;
       for (casadi_int k=0;k<n-1;++k) ind+=index_counters[k]*new_cumprod[k];

       // Fill in mapping
       for (casadi_int j=0;j<N_inner;++j) {
         mapping.at(m_ind++) = ind;
         ind+=stride;
       }

       // Bump first counter
       index_counters[0]++;

       // Overflow counters when needed
       for (casadi_int k=0;k<n-2;++k) {
         if (index_counters[k]==new_dims[k]) {
           index_counters[k] = 0;
           index_counters[k+1]++;
         }
       }
     }
     return mapping;
  }

  bvec_t* get_bvec_t(std::vector<double>& v) {
    if (v.empty()) {
      return nullptr;
    } else {
      return reinterpret_cast<bvec_t*>(&v.front());
    }
  }

  /// Get an pointer of sets of booleans from a double vector
  const bvec_t* get_bvec_t(const std::vector<double>& v) {
    if (v.empty()) {
      return nullptr;
    } else {
      return reinterpret_cast<const bvec_t*>(&v.front());
    }
  }

  std::string join(const std::vector<std::string>& l, const std::string& delim) {
    std::stringstream ss;
    for (casadi_int i=0;i<l.size();++i) {
      if (i>0) ss << delim;
      ss << l[i];
    }
    return ss.str();
  }

  bool startswith(const std::string& s, const std::string& p) {
    if (p.size()>s.size()) return false;
    for (casadi_int i=0;i<p.size();++i) {
      if (s[i]!=p[i]) return false;
    }
    return true;
  }

#ifdef HAVE_SIMPLE_MKSTEMPS
int simple_mkstemps_fd(const std::string& prefix, const std::string& suffix, std::string &result) {
    // Characters available for inventing filenames
    std::string chars = "abcdefghijklmnopqrstuvwxyz0123456789";
    int char_size = static_cast<int>(chars.size());

    // How many tries do we allow?
    casadi_int max_tries = std::numeric_limits<int>::max();

    // How long should the ID be to cover all tries?
    double max_tries_d = static_cast<double>(max_tries);
    double char_size_d = static_cast<double>(char_size);
    int id_size = lround(ceil(log(max_tries_d)/log(char_size_d)));

    // Random number generator
    std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<> r(0, char_size-1);

    for (casadi_int i=0;i<max_tries;++i) {
      result = prefix;
      for (casadi_int j=0;j<id_size;++j) {
        result += chars.at(r(rng));
      }
      result += suffix;

#ifdef _WIN32
      int fd = _sopen(result.c_str(),
        _O_BINARY | _O_CREAT | _O_EXCL | _O_RDWR, _SH_DENYNO, _S_IREAD | _S_IWRITE);
      // Could add _O_TEMPORARY, but then no possiblity of !cleanup_
#else
      int fd = open(result.c_str(), O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
#endif
      if (fd != -1) return fd;
      if (fd == -1 && errno != EEXIST) return -1;
    }
  }
std::string simple_mkstemps(const std::string& prefix, const std::string& suffix) {
  std::string ret;
  int fd = simple_mkstemps_fd(prefix, suffix, ret);
  if (fd==-1) {
    casadi_error("Failed to create temporary file: '" + ret + "'");
  } else {
#ifdef _WIN32
      _close(fd);
#else
      close(fd);
#endif
  }
  return ret;
}
#endif // HAVE_SIMPLE_MKSTEMPS

  std::string temporary_file(const std::string& prefix, const std::string& suffix) {
    #ifdef HAVE_MKSTEMPS
    // Preferred solution
    string ret = prefix + "XXXXXX" + suffix;
    if (mkstemps(&ret[0], static_cast<int>(suffix.size())) == -1) {
      casadi_error("Failed to create temporary file: '" + ret + "'");
    }
    return ret;
    #else // HAVE_MKSTEMPS
    #ifdef HAVE_SIMPLE_MKSTEMPS
    return simple_mkstemps(prefix, suffix);
    #else // HAVE_SIMPLE_MKSTEMPS
    // Fallback, may result in deprecation warnings
    return prefix + string(tmpnam(nullptr)) + suffix;
    #endif // HAVE_SIMPLE_MKSTEMPS
    #endif // HAVE_MKSTEMPS
  }


} // namespace casadi
