/*------------------------------------------------------------------------*/
/*! \file talisman.hpp
    \brief public library interface of TalisMan

  This is the only header needed to embed TalisMan into another application
  (e.g. as an application engine of Mallob). It deliberately does not include
  any internal header and only requires C++17.

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_INCLUDE_TALISMAN_TALISMAN_HPP_
#define TALISMAN_INCLUDE_TALISMAN_TALISMAN_HPP_
/*------------------------------------------------------------------------*/
#include <functional>
#include <memory>
#include <string>

namespace talisman {

/**
    Returns the version string of TalisMan, e.g. "2.0"
*/
const char*
version();

/*------------------------------------------------------------------------*/
// SAT backend
/*------------------------------------------------------------------------*/

/**
    Incremental SAT solver interface used for all SAT queries of TalisMan
    (IPASIR-like subset). Every solver instance is used by a single thread
    only, but several instances are used concurrently.
*/
class SatSolver {
  public:
  virtual ~SatSolver() = default;

  // Adds a literal to the current clause, 0 terminates the clause
  virtual void add(int lit) = 0;

  // Adds an assumption for the next call of solve()
  virtual void assume(int lit) = 0;

  // Returns 10 (SAT), 20 (UNSAT) or 0 (unknown, e.g. interrupted)
  virtual int solve() = 0;

  // After solve() returned 10: positive iff the literal is true
  virtual int val(int lit) = 0;
};

/**
    Creates a fresh solver instance. Called concurrently from worker threads,
    hence has to be thread-safe.
*/
using SatSolverFactory = std::function<std::unique_ptr<SatSolver>()>;

/**
    Creates the default (in-process CaDiCaL) solver. When used within run(),
    it honors the terminate callback given in Hooks.
*/
std::unique_ptr<SatSolver>
new_cadical_solver();

/*------------------------------------------------------------------------*/
// Running TalisMan
/*------------------------------------------------------------------------*/

enum class SpecMode { None, Miter, Mult, Assert, File };

struct Options {
  std::string input_file;// AIG to verify
  std::string spec_file; // only used with SpecMode::File
  SpecMode spec = SpecMode::None;
  int verbose = 1;          // 0..4
  int sc_depth = 2;         // subcircuit depth (levels)
  int primes = 0;           // 0: structured, 1: random
  int reduction = 0;        // 0: hybrid, 1: linear, 2: nonlinear
  bool custom_linalg = true;// custom routines for linear algebra
  bool uniform_sampling = true;
  bool trust_abc = false;// trust linear polynomial cache (skip SAT checks)
};

struct Hooks {
  // Factory for all SAT solvers, empty selects new_cadical_solver()
  SatSolverFactory sat_solver_factory;

  // Polled regularly (also from worker threads, hence has to be thread-safe);
  // returning true aborts the run with Status::Interrupted
  std::function<bool()> terminate;

  // Install TalisMan's own SIGINT/SIGSEGV/SIGABRT/SIGTERM handlers. Only
  // useful for the command line front end, embedding applications keep
  // their own handlers.
  bool install_signal_handlers = false;
};

enum class Status {
  Verified,   // specification holds
  NotVerified,// guess and prove failed
  Interrupted,// terminate callback fired
  Error       // invalid input or internal error, see error_code
};

struct Result {
  Status status = Status::Error;
  int error_code = 0;// error code of the failing check, 255 if unexpected
  std::string error_message;
};

/**
    Verifies the given AIG against the selected specification.

    TalisMan keeps its state in global variables, hence concurrent calls are
    serialized. All state is reset between runs, so run() can be called
    repeatedly within the same process.
*/
Result
run(const Options& options, const Hooks& hooks = Hooks());

}// namespace talisman

#endif// TALISMAN_INCLUDE_TALISMAN_TALISMAN_HPP_
