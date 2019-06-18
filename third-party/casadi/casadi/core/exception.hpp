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


#ifndef CASADI_EXCEPTION_HPP
#define CASADI_EXCEPTION_HPP

#include <exception>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <iomanip>
#include <chrono>

#include <casadi/core/casadi_export.h>

namespace casadi {

/** \brief  Casadi exception class
        \author Joel Andersson
        \date 2010
        Example for simple exception throwing:
        \code
                throw CasadiException("This is a nasty error");
        \endcode
        Example for exception chaining:
        \code
                try {
                        throw CasadiException("This is a nasty error");
                catch(CasadiException &e) {
                        throw CasadiException("Serious error.") << e;
                }
        \endcode
*/
class CasadiException : public std::exception {
  public:
  //! \brief Default constructor
  CasadiException() {
  }

  //! \brief Form message string
  explicit CasadiException(const std::string& msg) : msg_(msg) {}

  //! \brief Destructor
  ~CasadiException() throw() {}

  //! \brief Display error
  const char* what() const throw() override {
    return msg_.c_str();
  }

  protected:
  std::string msg_;
};

class KeyboardInterruptException : public CasadiException {
  public:
  //! \brief Default constructor
  KeyboardInterruptException() : CasadiException("KeyboardInterrupt") {}
  //! \brief Destructor
  ~KeyboardInterruptException() throw() {}
};

// Strip path prefix
inline std::string trim_path(const std::string& full_path) {
  size_t found = full_path.rfind("/casadi/");
  if (found == std::string::npos) {
    return full_path;
  } else {
    std::string ret = full_path;
    ret.replace(0, found, "...");
    return ret;
  }
}

// Current time as a string
inline std::ostream& message_prefix(std::ostream &stream) {
  stream << "CasADi - ";
  auto time = std::time(nullptr);
  char stamp[30];
  strftime(stamp, 30, "%F %T", std::localtime(&time)); // NOLINT(runtime/threadsafe_fn)
  stream << stamp;
  return stream;
}

// String denoting where the macro is situated
#define CASADI_WHERE casadi::trim_path(__FILE__ ":" CASADI_STR(__LINE__))

// Throw an exception with information about source code location
#define casadi_error(msg, ...) \
throw casadi::CasadiException(CASADI_WHERE + ": "\
          + casadi::fmtstr(msg, casadi::strvec(__VA_ARGS__)))

// This assertion checks for illegal user inputs
#define casadi_assert(x, msg, ...) \
if (!(x)) casadi_error("Assertion \"" CASADI_STR(x) "\" failed:\n"\
          + std::string(msg), __VA_ARGS__)

// This assertion if for internal errors caused by bugs in CasADi
#define casadi_assert_dev(x) casadi_assert(x, "Notify the CasADi developers.")

// This assertion if for internal errors caused by bugs in CasADi
#define casadi_report() casadi_error("Notify the CasADi developers.")

// Issue a warning, including location in the source code
#define casadi_warning(msg) \
  casadi::message_prefix(casadi::uerr()) \
    << " WARNING(\"" << msg << "\") [" << CASADI_WHERE << "]\n" << std::flush;

// Issue a message, including location in the source code
#define casadi_message(msg) \
  casadi::message_prefix(casadi::uout()) \
    << " MESSAGE(\"" << msg << "\") [" << CASADI_WHERE << "]\n" << std::flush;

} // namespace casadi

#endif // CASADI_EXCEPTION_HPP
