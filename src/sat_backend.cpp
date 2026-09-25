/*------------------------------------------------------------------------*/
/*! \file sat_backend.cpp
    \brief selects the SAT solver used for all SAT queries

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "sat_backend.hpp"

#include <cadical.hpp>

#include "run_control.hpp"

namespace talisman {

static SatSolverFactory sat_solver_factory;

/*------------------------------------------------------------------------*/
// Default backend: in-process CaDiCaL
class CadicalSolver final : public SatSolver, private CaDiCaL::Terminator {
  public:
  CadicalSolver() {
    // only pay for the terminator checks if someone can actually terminate
    if(has_terminate_callback())
      solver.connect_terminator(this);
  }
  ~CadicalSolver() override { solver.disconnect_terminator(); }

  void add(int lit) override { solver.add(lit); }
  void assume(int lit) override { solver.assume(lit); }
  int solve() override { return solver.solve(); }
  int val(int lit) override { return solver.val(lit); }

  private:
  bool terminate() override { return terminate_requested(); }

  CaDiCaL::Solver solver;
};

/*------------------------------------------------------------------------*/
std::unique_ptr<SatSolver>
new_cadical_solver() {
  return std::make_unique<CadicalSolver>();
}
/*------------------------------------------------------------------------*/
void
set_sat_solver_factory(SatSolverFactory factory) {
  sat_solver_factory = std::move(factory);
}
/*------------------------------------------------------------------------*/
std::unique_ptr<SatSolver>
new_sat_solver() {
  if(sat_solver_factory)
    return sat_solver_factory();
  return new_cadical_solver();
}

}// namespace talisman
