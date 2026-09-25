/*------------------------------------------------------------------------*/
/*! \file sat_backend.hpp
    \brief selects the SAT solver used for all SAT queries

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SAT_BACKEND_H_
#define TALISMAN_SRC_SAT_BACKEND_H_
/*------------------------------------------------------------------------*/
#include "talisman/talisman.hpp"

namespace talisman {

/**
    Installs the factory used by new_sat_solver(), empty selects CaDiCaL.
    Must not be called while solvers are created.
*/
void
set_sat_solver_factory(SatSolverFactory factory);

/**
    Creates a fresh solver from the installed factory. Thread-safe.
*/
std::unique_ptr<SatSolver>
new_sat_solver();

}// namespace talisman

#endif// TALISMAN_SRC_SAT_BACKEND_H_
