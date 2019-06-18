cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Conic nlpsol
  qp_to_nlp.hpp qp_to_nlp.cpp qp_to_nlp_meta.cpp)

# Active-set QP solver
casadi_plugin(Conic qrqp qrqp.hpp qrqp.cpp qrqp_meta.cpp)

# Simple just-in-time compiler, using shell commands
if(WITH_DL)
  casadi_plugin(Importer shell
    shell_compiler.hpp
    shell_compiler.cpp
    shell_compiler_meta.cpp)
endif()

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  runge_kutta.hpp
  runge_kutta.cpp
  runge_kutta_meta.cpp)

# Collocation integrator
casadi_plugin(Integrator collocation
  collocation.hpp
  collocation.cpp
  collocation_meta.cpp)

# Linear interpolant
casadi_plugin(Interpolant linear
  linear_interpolant.hpp linear_interpolant.cpp linear_interpolant_meta.cpp
)

casadi_plugin(Interpolant bspline
  bspline_interpolant.hpp bspline_interpolant.cpp bspline_interpolant_meta.cpp
)

casadi_plugin(Linsol symbolicqr
  symbolic_qr.hpp symbolic_qr.cpp symbolic_qr_meta.cpp
)

# Sparse direct QR - implemented in CasADi's C runtime
casadi_plugin(Linsol qr
  linsol_qr.hpp linsol_qr.cpp linsol_qr_meta.cpp
)

# Sparse direct LDL' - implemented in CasADi's C runtime
casadi_plugin(Linsol ldl
  linsol_ldl.hpp linsol_ldl.cpp linsol_ldl_meta.cpp
)

casadi_plugin(Linsol lsqr
  lsqr.hpp lsqr.cpp lsqr_meta.cpp
)

# SQPMethod -  A basic SQP method
casadi_plugin(Nlpsol sqpmethod
  sqpmethod.hpp sqpmethod.cpp sqpmethod_meta.cpp)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_plugin(Nlpsol scpgen
  scpgen.hpp scpgen.cpp scpgen_meta.cpp)

casadi_plugin(Rootfinder newton
  newton.hpp newton.cpp newton_meta.cpp)

casadi_plugin(Rootfinder fast_newton
  fast_newton.hpp fast_newton.cpp fast_newton_meta.cpp)

casadi_plugin(Rootfinder nlpsol
  implicit_to_nlp.hpp implicit_to_nlp.cpp implicit_to_nlp_meta.cpp)
