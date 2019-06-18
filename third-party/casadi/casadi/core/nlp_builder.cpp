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


#include "nlp_builder.hpp"
#include "core.hpp"
#include <fstream>

using namespace std;
namespace casadi {

  void NlpBuilder::import_nl(const std::string& filename, const Dict& opts) {
    // Redirect to helper class
    NlImporter(*this, filename, opts);
  }

  void NlpBuilder::disp(std::ostream& stream, bool more) const {
    stream << "#x=" << this->x.size() << ", #g=" << this->g.size();
    if (more) {
      stream << endl;
      stream << "x = " << this->x << endl;
      stream << "f = " << this->f << endl;
      stream << "g = " << this->g << endl;
    }
  }

  NlImporter::NlImporter(NlpBuilder& nlp, const std::string& filename, const Dict& opts)
  : nlp_(nlp) {
    // Set default options
    verbose_=false;

    // Read user options
    for (auto&& op : opts) {
      if (op.first == "verbose") {
        verbose_ = op.second;
      } else {
        stringstream ss;
        ss << "Unknown option \"" << op.first << "\"" << endl;
        throw CasadiException(ss.str());
      }
    }
    // Open file for reading
    s_.open(filename.c_str());
    if (verbose_) casadi_message("Reading file \"" + filename + "\"");

    // Read the header of the NL-file (first 10 lines)
    const casadi_int header_sz = 10;
    vector<string> header(header_sz);
    for (casadi_int k=0; k<header_sz; ++k) {
      getline(s_, header[k]);
    }

    // Assert that the file is not in binary form
    if (header.at(0).at(0)=='g') {
      binary_ = false;
    } else if (header.at(0).at(0)=='b') {
      binary_ = true;
    } else {
      casadi_error("File could not be read");
    }

    // Get the number of objectives and constraints
    stringstream ss(header[1]);
    ss >> n_var_ >> n_con_ >> n_obj_ >> n_eq_ >> n_lcon_;
    if (verbose_) {
      casadi_message("n_var=" + str(n_var_) + ", n_con =" + str(n_con_) + ", "
                     "n_obj=" + str(n_obj_) + ", n_eq=" + str(n_eq_) + ", "
                     "n_lcon=" + str(n_lcon_));
    }

    // Get the number of nonlinear vars in constraints, objectives, both
    stringstream ss4(header[4]);
    ss4 >> nlvc_ >> nlvo_ >> nlvb_;
    if (verbose_) {
      casadi_message("nlvc=" + str(nlvc_) + ", nlvo=" + str(nlvo_) + ", nlvb=" + str(nlvb_));
    }

    // Get the number of discrete variables
    stringstream ss6(header[6]);
    ss6 >> nbv_ >> niv_ >> nlvbi_ >> nlvci_ >> nlvoi_;
    if (verbose_) {
      casadi_message("nbv=" + str(nbv_) + ", niv =" + str(niv_) + ", "
                     "nlvbi=" + str(nlvbi_) + ", nlvci=" + str(nlvci_) + ", "
                     "nlvoi=" + str(nlvoi_));
    }

    // Allocate variables
    nlp_.x = MX::sym("x", 1, 1, n_var_);

    // Allocate f and c
    nlp_.f = 0;
    nlp_.g.resize(n_con_, 0);

    // Allocate bounds for x and primal initial guess
    nlp_.x_lb.resize(n_var_, -inf);
    nlp_.x_ub.resize(n_var_,  inf);
    nlp_.x_init.resize(n_var_, 0);

    // Allocate bounds for g and dual initial guess
    nlp_.g_lb.resize(n_con_, -inf);
    nlp_.g_ub.resize(n_con_,  inf);
    nlp_.lambda_init.resize(n_con_, 0);

    // Allocate binary variables vector
    nlp_.discrete.clear();

    //D. M. Gay and M. Hill, 'Hooking Your Solver to AMPL' October, 1997.
    // continuous in an objective and in a constraint
    for (casadi_int j=0; j<nlvb_-nlvbi_; ++j) nlp_.discrete.push_back(false);

    // integer in an objective and in a constraint
    for (casadi_int j=0; j<nlvbi_; ++j) nlp_.discrete.push_back(true);

    // continuous just in constraints
    for (casadi_int j=0; j<nlvc_ - (nlvb_ + nlvci_); ++j) nlp_.discrete.push_back(false);

    // integer just in constraints
    for (casadi_int j=0; j<nlvci_; ++j) nlp_.discrete.push_back(true);

    // continuous just in objectives
    for (casadi_int j=0; j<nlvo_ - (nlvc_ + nlvoi_); ++j) nlp_.discrete.push_back(false);

    // integer just in objectives
    for (casadi_int j=0; j < nlvoi_; ++j) nlp_.discrete.push_back(true);

    // linear
    casadi_int max_nlvc_nlvo = (nlvc_ < nlvo_) ? nlvo_ : nlvc_;
    for (casadi_int j=0; j<n_var_-(max_nlvc_nlvo+niv_+nbv_); ++j) nlp_.discrete.push_back(false);

    // binary
    for (casadi_int j = 0; j<nbv_; ++j) nlp_.discrete.push_back(true);

    // other integer
    for (casadi_int j = 0; j<niv_; ++j) nlp_.discrete.push_back(true);

    casadi_assert(nlp_.discrete.size()==n_var_,
      "Number of variables in the header don't match");

    // All variables, including dependent
    v_ = nlp_.x;

    if (binary_) {
      streampos offset = s_.tellg();
      s_.close();
      s_.open(filename.c_str(), std::ifstream::binary);
      s_.seekg(offset);
    }

    // Read segments
    parse();

    // multiple the objective sign
    nlp_.f = sign_*nlp_.f;
  }

  NlImporter::~NlImporter() {
    // Close the NL file
    s_.close();
  }

  void NlImporter::parse() {
    // Segment key
    char key;

    // Process segments
    while (true) {
      // Read segment key
      key = read_char();
      if (s_.eof()) break; // end of file encountered
      switch (key) {
        case 'F': F_segment(); break;
        case 'S': S_segment(); break;
        case 'V': V_segment(); break;
        case 'C': C_segment(); break;
        case 'L': L_segment(); break;
        case 'O': O_segment(); break;
        case 'd': d_segment(); break;
        case 'x': x_segment(); break;
        case 'r': r_segment(); break;
        case 'b': b_segment(); break;
        case 'k': k_segment(); break;
        case 'J': J_segment(); break;
        case 'G': G_segment(); break;
        default: casadi_error("Unknown .nl segment");
      }
    }
  }

  MX NlImporter::expr() {
    // Read the instruction
    char inst = read_char();

    // Temporaries
    int i;
    double d;

    // Error message
    stringstream msg;

    // Process instruction
    switch (inst) {

      // Symbolic variable
      case 'v':
      // Read the variable number
      i = read_int();

      // Return the corresponding expression
      return v_.at(i);

      // Numeric expression
      case 'n':

      // Read the floating point number
      d = read_double();

      // Return an expression containing the number
      return d;

      // Numeric expression
      case 's':

      // Read the short number
      d = read_short();

      // Return an expression containing the number
      return d;

      // Numeric expression
      case 'l':

      // Read the short number
      d = static_cast<double>(read_long());

      // Return an expression containing the number
      return d;

      // Operation
      case 'o':

      // Read the operation
      i = read_int();

      // Process
      switch (i) {

        // Unary operations, class 1 in Gay2005
        case 13:  case 14:  case 15:  case 16:  case 34:  case 37:  case 38:  case 39:  case 40:
        case 41:  case 43:  case 42:  case 44:  case 45:  case 46:  case 47:  case 49:  case 50:
        case 51:  case 52:  case 53:
        {
          // Read dependency
          MX x = expr();

          // Perform operation
          switch (i) {
            case 13:  return floor(x);
            case 14:  return ceil(x);
            case 15:  return abs(x);
            case 16:  return -x;
            case 34:  return logic_not(x);
            case 37:  return tanh(x);
            case 38:  return tan(x);
            case 39:  return sqrt(x);
            case 40:  return sinh(x);
            case 41:  return sin(x);
            case 42:  return log10(x);
            case 43:  return log(x);
            case 44:  return exp(x);
            case 45:  return cosh(x);
            case 46:  return cos(x);
            // case 47:  return atanh(x); FIXME
            case 49:  return atan(x);
            // case 50:  return asinh(x); FIXME
            case 51:  return asin(x);
            // case 52:  return acosh(x); FIXME
            case 53:  return acos(x);

            default:
            casadi_error("Unknown unary operation: " + str(i));
          }
          break;
        }

        // Binary operations, class 2 in Gay2005
        case 0:   case 1:   case 2:   case 3:   case 4:   case 5:   case 6:   case 20:  case 21:
        case 22:  case 23:  case 24:  case 28:  case 29:  case 30:  case 48:  case 55:  case 56:
        case 57:  case 58:  case 73:
        {
          // Read dependencies
          MX x = expr();
          MX y = expr();

          // Perform operation
          switch (i) {
            case 0:   return x + y;
            case 1:   return x - y;
            case 2:   return x * y;
            case 3:   return x / y;
            // case 4:   return rem(x, y); FIXME
            case 5:   return pow(x, y);
            // case 6:   return x < y; // TODO(Joel): Verify this,
            // what is the difference to 'le' == 23 below?
            case 20:  return logic_or(x, y);
            case 21:  return logic_and(x, y);
            case 22:  return x < y;
            case 23:  return x <= y;
            case 24:  return x == y;
            case 28:  return x >= y;
            case 29:  return x > y;
            case 30:  return x != y;
            case 48:  return atan2(x, y);
            // case 55:  return intdiv(x, y); // FIXME
            // case 56:  return precision(x, y); // FIXME
            // case 57:  return round(x, y); // FIXME
            // case 58:  return trunc(x, y); // FIXME
            // case 73:  return iff(x, y); // FIXME

            default:
            casadi_error("Unknown binary operation: " + str(i));
          }
          break;
        }

        // N-ary operator, classes 2, 6 and 11 in Gay2005
        case 11: case 12: case 54: case 59: case 60: case 61: case 70: case 71: case 74:
        {
          // Number of elements in the sum
          int n = read_int();

          // Collect the arguments
          vector<MX> args(n);
          for (int k=0; k<n; ++k) {
            args[k] = expr();
          }

          // Perform the operation
          switch (i) {
            // case 11: return min(args).scalar(); FIXME // rename?
            // case 12: return max(args).scalar(); FIXME // rename?
            // case 54: return sum(args).scalar(); FIXME // rename?
            // case 59: return count(args).scalar(); FIXME // rename?
            // case 60: return numberof(args).scalar(); FIXME // rename?
            // case 61: return numberofs(args).scalar(); FIXME // rename?
            // case 70: return all(args).scalar(); FIXME // and in AMPL // rename?
            // case 71: return any(args).scalar(); FIXME // or in AMPL // rename?
            // case 74: return alldiff(args).scalar(); FIXME // rename?
            case 54:
            {
              MX r = 0;
              for (vector<MX>::const_iterator it=args.begin();
              it!=args.end(); ++it) r += *it;
              return r;
            }

            default:
            casadi_error("Unknown n-ary operation: " + str(i));
          }
          break;
        }

        // Piecewise linear terms, class 4 in Gay2005
        case 64:
        casadi_error("Piecewise linear terms not supported");
        break;

        // If-then-else expressions, class 5 in Gay2005
        case 35: case 65: case 72:
        casadi_error("If-then-else expressions not supported");
        break;

        default:
        casadi_error("Unknown operation: " + str(i));
      }
      break;

      default:
       uout() << s_.tellg() << std::endl;
      casadi_error("Unknown instruction: " + str(inst));
    }

    // Throw error message
    casadi_error("Unknown error");
  }

  void NlImporter::F_segment() {
    casadi_error("Imported function description unsupported.");
  }

  void NlImporter::S_segment() {
    casadi_error("Suffix values unsupported");
  }

  void NlImporter::V_segment() {
    // Read header
    int i = read_int();
    int j = read_int();
    read_int();

    // Make sure that v is long enough
    if (i >= v_.size()) {
      v_.resize(i+1);
    }

    // Initialize element to zero
    v_.at(i) = 0;

    // Add the linear terms
    for (int jj=0; jj<j; ++jj) {
      // Linear term
      int pl = read_int();
      double cl = read_double();

      // Add to variable definition (assuming it has already been defined)
      casadi_assert(!v_.at(pl).is_empty(), "Circular dependencies not supported");
      v_.at(i) += cl*v_.at(pl);
    }

    // Finally, add the nonlinear term
    v_.at(i) += expr();
  }

  int NlImporter::read_int() {
    int i;
    if (binary_) {
      s_.read(reinterpret_cast<char *>(&i), sizeof(int));
    } else {
      s_ >> i;
    }
    return i;
  }

  char NlImporter::read_char() {
    char c;
    if (binary_) {
      s_.read(&c, 1);
    } else {
      s_ >> c;
    }
    return c;
  }

  double NlImporter::read_double() {
    double d;
    if (binary_) {
      s_.read(reinterpret_cast<char *>(&d), sizeof(double));
    } else {
      s_ >> d;
    }
    return d;
  }

  short NlImporter::read_short() {
    short d;
    if (binary_) {
      s_.read(reinterpret_cast<char *>(&d), 2);
    } else {
      s_ >> d;
    }
    return d;
  }

  long NlImporter::read_long() {
    long d;
    if (binary_) {
      s_.read(reinterpret_cast<char *>(&d), 4);
    } else {
      s_ >> d;
    }
    return d;
  }

  void NlImporter::C_segment() {
    // Get the number
    int i = read_int();

    // Parse and save expression
    nlp_.g.at(i) = expr();
  }

  void NlImporter::L_segment() {
    casadi_error("Logical constraint expression unsupported");
  }

  void NlImporter::O_segment() {
    // Get the number
    read_int(); // i

    // Should the objective be maximized
    int sigma= read_int();
    sign_ = sigma!=0 ? -1 : 1;

    // Parse and save expression
    nlp_.f += expr();
  }

  void NlImporter::d_segment() {
    // Read the number of guesses supplied
    int m = read_int();

    // Process initial guess for the fual variables
    for (int i=0; i<m; ++i) {
      // Offset and value
      int offset = read_int();
      double d = read_double();

      // Save initial guess
      nlp_.lambda_init.at(offset) = d;
    }
  }

  void NlImporter::x_segment() {
    // Read the number of guesses supplied
    int m = read_int();

    // Process initial guess
    for (int i=0; i<m; ++i) {
      // Offset and value
      int offset = read_int();
      double d = read_double();

      // Save initial guess
      nlp_.x_init.at(offset) = d;
    }
  }

  void NlImporter::r_segment() {
    // For all constraints
    for (int i=0; i<n_con_; ++i) {

      // Read constraint type
      char c_type = read_char();

      // Temporary
      double c;

      switch (c_type) {
        // Upper and lower bounds
        case '0':
        c = read_double();
        nlp_.g_lb.at(i) = c;
        c = read_double();
        nlp_.g_ub.at(i) = c;
        continue;

        // Only upper bounds
        case '1':
        c = read_double();
        nlp_.g_ub.at(i) = c;
        continue;

        // Only lower bounds
        case '2':
        c = read_double();
        nlp_.g_lb.at(i) = c;
        continue;

        // No bounds
        case '3':
        continue;

        // Equality constraints
        case '4':
        c = read_double();
        nlp_.g_lb.at(i) = nlp_.g_ub.at(i) = c;
        continue;

        // Complementary constraints
        case '5':
        {
          // Read the indices
          read_int(); // ck
          read_int(); // ci
          casadi_error("Complementary constraints unsupported");
          continue;
        }

        default:
        casadi_error("Illegal constraint type");
      }
    }
  }

  void NlImporter::b_segment() {
    // For all variable
    for (casadi_int i=0; i<n_var_; ++i) {

      // Read constraint type
      char c_type = read_char();

      // Temporary
      double c;

      switch (c_type) {
        // Upper and lower bounds
        case '0':
        c = read_double();
        nlp_.x_lb.at(i) = c;
        c = read_double();
        nlp_.x_ub.at(i) = c;
        continue;

        // Only upper bounds
        case '1':
        c = read_double();
        nlp_.x_ub.at(i) = c;
        continue;

        // Only lower bounds
        case '2':
        c = read_double();
        nlp_.x_lb.at(i) = c;
        continue;

        // No bounds
        case '3':
        continue;

        // Equality constraints
        case '4':
        c = read_double();
        nlp_.x_lb.at(i) = nlp_.x_ub.at(i) = c;
        continue;

        default:
        casadi_error("Illegal variable bound type");
      }
    }
  }

  void NlImporter::k_segment() {
    // Get row offsets
    vector<casadi_int> rowind(n_var_+1);

    // Get the number of offsets
    int k = read_int();
    casadi_assert_dev(k==n_var_-1);

    // Get the row offsets
    rowind[0]=0;
    for (int i=0; i<k; ++i) {
      rowind[i+1] = read_int();
    }
  }

  void NlImporter::J_segment() {
    // Get constraint number and number of terms
    int i = read_int();
    int k = read_int();

    // Get terms
    for (int kk=0; kk<k; ++kk) {
      // Get the term
      int j = read_int();
      double c = read_double();

      // Add to constraints
      nlp_.g.at(i) += c*v_.at(j);
    }
  }

  void NlImporter::G_segment() {
    // Get objective number and number of terms
    read_int(); // i
    int k = read_int();

    // Get terms
    for (int kk=0; kk<k; ++kk) {
      // Get the term
      int j = read_int();
      double c = read_double();

      // Add to objective
      nlp_.f += c*v_.at(j);
    }
  }


} // namespace casadi
