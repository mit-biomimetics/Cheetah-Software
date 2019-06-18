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


#include "options.hpp"
#include <algorithm>
#include <locale>

using namespace std;

namespace casadi {

  const Options::Entry* Options::find(const std::string& name) const {
    // Check if in one of the bases
    for (auto&& b : bases) {
      const Options::Entry* entry = b->find(name);
      if (entry) return entry;
    }

    // Lookup in this class
    auto it = entries.find(name);
    if (it!=entries.end()) {
      return &it->second;
    } else {
      return nullptr;
    }
  }

  void Options::Entry::disp(const std::string& name, std::ostream &stream) const {
    stream << "> \"" << name << "\"          ["
           << GenericType::get_type_description(this->type)
           << "] ";

    // Print out the description on a new line.
    stream << "     \"" << this->description << "\""<< std::endl;
  }

  void Options::disp(std::ostream& stream) const {
    // Print bases
    for (auto&& b : bases) {
      b->disp(stream);
    }

    // Print all entries
    for (auto&& e : entries) {
      e.second.disp(e.first, stream);
    }
  }

  double Options::word_distance(const std::string &a, const std::string &b) {
    /// Levenshtein edit distance
    if (a == b) return 0;
    casadi_int na = a.size();
    casadi_int nb = b.size();
    if (na == 0) return static_cast<double>(nb);
    if (nb == 0) return static_cast<double>(na);

    vector<casadi_int> v0(nb+1, 0);
    vector<casadi_int> v1(nb+1, 0);

    for (casadi_int i=0;i<nb+1;++i)
      v0[i] = i;

    char s;
    char t;
    std::locale loc;
    for (casadi_int i=0;i<na;i++) {
      v1[0] = i + 1;
      for (casadi_int j=0; j<nb; j++) {
        s = std::tolower(a[i], loc);
        t = std::tolower(b[j], loc);
        casadi_int cost = 0;
        if (s != t)
          cost = 1;

        v1[j+1] = min(min(v1[j] + 1, v0[j+1] + 1), v0[j] + cost);
      }

      for (casadi_int j=0; j<nb+1; j++)
        v0[j] = v1[j];
    }

    return static_cast<double>(v1[nb]);
  }

  vector<string> Options::suggestions(const string& word, casadi_int amount) const {
    // Best distances so far
    const double inf = numeric_limits<double>::infinity();
    vector<pair<double, string> > best(amount, {inf, ""});

    // Iterate over elements
    best_matches(word, best);

    // Sort the elements in ascending order
    stable_sort(best.begin(), best.end());

    // Collect the values that are non-infinite
    vector<string> ret;
    ret.reserve(amount);
    for (auto&& e : best) {
      if (e.first!=inf) {
        ret.push_back(e.second);
      }
    }
    return ret;
  }

  void Options::best_matches(const std::string& word,
                             vector<pair<double, string> >& best) const {
    // Iterate over bases
    for (auto&& b : bases) {
      b->best_matches(word, best);
    }

    // Worst match so far
    auto worst = max_element(best.begin(), best.end());

    // Loop over entries
    for (auto&& e : entries) {
      // Get word distance
      double d = word_distance(e.first, word);

      // Keep if better than the worst amongst the suggestions
      if (d < worst->first) {
        worst->first = d;
        worst->second = e.first;
        worst = max_element(best.begin(), best.end());
      }
    }
  }

  bool Options::has_dot(const Dict& opts) {
    for (auto&& op : opts) {
      if (op.first.find('.') != string::npos || op.first.find("__") != string::npos) {
        return true;
      }
    }
    return false;
  }

  bool Options::has_null(const Dict& opts) {
    for (auto&& op : opts) {
      if (op.second.is_null()) return true;
    }
    return false;
  }

  bool Options::is_sane(const Dict& opts) {
    return !has_dot(opts) && !has_null(opts);
  }

  Dict Options::sanitize(const Dict& opts) {
    // Drop nulls
    if (has_null(opts)) {
      // Create a new dictionary without the null entries
      Dict ret;
      for (auto&& op : opts) {
        if (!op.second.is_null()) ret[op.first] = op.second;
      }
      return ret;
    }

    //  Treat the case where any of the options have a dot (dictionary shorthand)
    if (has_dot(opts)) {
      // New options dictionary being constructed
      Dict ret;

      // Sub-dictionary and corresponding name being constructed
      Dict sopts;
      std::string sname;

      // Process options
      for (auto&& op : opts) {
        // Find the dot if any
        string::size_type dotpos = op.first.find('.'), dotpos_end;
        if (dotpos==string::npos) {
          dotpos = op.first.find("__");
          if (dotpos!=string::npos) dotpos_end = dotpos+2;
        } else {
          dotpos_end = dotpos+1;
        }

        // Flush last sub-dictionary
        if (!sname.empty() && (dotpos==string::npos
                               || op.first.compare(0, dotpos, sname)!=0)) {
          ret[sname] = sopts;
          sname.clear();
          sopts.clear();
        }

        // Add to dictionary
        if (dotpos != string::npos) {
          sname = op.first.substr(0, dotpos);
          sopts[op.first.substr(dotpos_end)] = op.second;
        } else {
          ret[op.first] = op.second;
        }
      }

      // Flush trailing sub-dictionary
      if (!sname.empty()) ret[sname] = sopts;

      return ret;
    }

    // Nothing to do
    return opts;
  }

  void Options::check(const Dict& opts) const {
    // Make sure all options exist and have the correct type
    for (auto&& op : opts) {
      const Options::Entry* entry = find(op.first);

      // Informative error message if option does not exist
      if (entry==nullptr) {
        stringstream ss;
        ss << "Unknown option: " << op.first << endl;
        ss << endl;
        ss << "Did you mean one of the following?" << endl;
        for (auto&& s : suggestions(op.first)) {
          print_one(s, ss);
        }
        ss << "Use print_options() to get a full list of options." << endl;
        casadi_error(ss.str());
      }

      // Check type
      casadi_assert(op.second.can_cast_to(entry->type),
                            "Illegal type for " + op.first + ": " +
                            op.second.get_description() +
                            " cannot be cast to " +
                            GenericType::get_type_description(entry->type) + ".");

    }
  }

  void Options::print_all(std::ostream &stream) const {
    stream << "\"Option name\" [type] = value" << endl;
    disp(stream);
    stream << endl;
  }

  void Options::print_one(const std::string &name, std::ostream &stream) const {
    const Options::Entry* entry = find(name);
    if (entry!=nullptr) {
      entry->disp(name, stream);
    } else {
      stream << "  \"" << name << "\" does not exist.";
    }
  }

  std::vector<std::string> Options::all() const {
    std::vector<std::string> ret;
    for (auto&& e : entries) ret.push_back(e.first);
    return ret;
  }

  std::string Options::type(const std::string& name) const {
    const Options::Entry* entry = find(name);
    casadi_assert(entry!=nullptr, "Option \"" + name + "\" does not exist");
    return GenericType::get_type_description(entry->type);
  }

  std::string Options::info(const std::string& name) const {
    const Options::Entry* entry = find(name);
    casadi_assert(entry!=nullptr, "Option \"" + name + "\" does not exist");
    return entry->description;
  }

} // namespace casadi
