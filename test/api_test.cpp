/*------------------------------------------------------------------------*/
/*! \file api_test.cpp
    \brief tests the library interface as used by embedding applications:
    repeated runs within one process, custom SAT backends, interruption and
    recovery after errors

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include <atomic>
#include <cstdio>
#include <memory>
#include <string>

#include "talisman/talisman.hpp"

using namespace talisman;

static int failures = 0;

#define CHECK(cond)                                              \
  do {                                                           \
    if(!(cond)) {                                                \
      std::fprintf(stderr, "FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      failures++;                                                \
    }                                                            \
  } while(0)

static Options
mult_options(const char* file) {
  Options opts;
  opts.input_file = std::string(TALISMAN_EXAMPLES_DIR) + "/" + file;
  opts.spec = SpecMode::Mult;
  opts.verbose = 0;
  return opts;
}

/*------------------------------------------------------------------------*/
// Forwards to CaDiCaL and counts the calls, optionally reporting "unknown"
static std::atomic<int> num_solvers{ 0 };
static std::atomic<int> num_solves{ 0 };

class CountingSolver : public SatSolver {
  public:
  explicit CountingSolver(bool answer_unknown)
    : inner(new_cadical_solver())
    , answer_unknown(answer_unknown) {
    num_solvers++;
  }
  void add(int lit) override { inner->add(lit); }
  void assume(int lit) override { inner->assume(lit); }
  int solve() override {
    num_solves++;
    return answer_unknown ? 0 : inner->solve();
  }
  int val(int lit) override { return inner->val(lit); }

  private:
  std::unique_ptr<SatSolver> inner;
  bool answer_unknown;
};

/*------------------------------------------------------------------------*/
int
main() {
  // repeated runs on different circuits within the same process
  CHECK(run(mult_options("bp-ar-cl-6.aig")).status == Status::Verified);
  CHECK(run(mult_options("abc8.aig")).status == Status::Verified);
  CHECK(run(mult_options("bp-ar-cl-6.aig")).status == Status::Verified);

  // all SAT queries go through the installed backend
  {
    Hooks hooks;
    hooks.sat_solver_factory = [] {
      return std::make_unique<CountingSolver>(false);
    };
    CHECK(run(mult_options("bp-ar-cl-6.aig"), hooks).status
          == Status::Verified);
    CHECK(num_solvers > 0);
    CHECK(num_solves > 0);
    std::printf("custom backend: %d solvers, %d solve calls\n",
                num_solvers.load(),
                num_solves.load());
  }

  // an unknown SAT answer must never be taken as a proof
  {
    Hooks hooks;
    hooks.sat_solver_factory = [] {
      return std::make_unique<CountingSolver>(true);
    };
    Result res = run(mult_options("bp-ar-cl-6.aig"), hooks);
    CHECK(res.status == Status::Error);
    CHECK(res.error_code == 8);
  }

  // interruption right away and in the middle of the run
  for(int polls : { 0, 200 }) {
    std::atomic<int> count{ 0 };
    Hooks hooks;
    hooks.terminate = [&count, polls] { return count++ >= polls; };
    CHECK(run(mult_options("bp-ar-cl-6.aig"), hooks).status
          == Status::Interrupted);
  }

  // errors are reported instead of exiting the process
  {
    Result res = run(mult_options("does-not-exist.aig"));
    CHECK(res.status == Status::Error);
    CHECK(res.error_code == 20);
  }

  // the process is still fully usable afterwards
  CHECK(run(mult_options("bp-ar-cl-6.aig")).status == Status::Verified);

  if(failures)
    std::fprintf(stderr, "%d check(s) failed\n", failures);
  else
    std::printf("all checks passed\n");
  return failures ? 1 : 0;
}
