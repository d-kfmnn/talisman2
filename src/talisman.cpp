/*------------------------------------------------------------------------*/
/*! \file talisman.cpp
    \brief library entry point of TalisMan, see include/talisman/talisman.hpp

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#define VERSION "2.0"
/*------------------------------------------------------------------------*/
#include <mutex>
#include <optional>

#include "talisman/talisman.hpp"

#include "candidates.hpp"
#include "coeff_vec.hpp"
#include "gate.hpp"
#include "guess_and_prove.hpp"
#include "linear_algebra.hpp"
#include "linpoly_cache.hpp"
#include "linpolynomial.hpp"
#include "parser.hpp"
#include "run_control.hpp"
#include "sat_backend.hpp"
#include "signal_statistics.hpp"
#include "specpoly.hpp"
#include "subcircuit.hpp"
#include "thread_pool.hpp"

using namespace talisman;

/*------------------------------------------------------------------------*/

static const int prime = 101;

// Define linpolynomial pruning flag (default enabled)
namespace talisman {
bool custom_linalg = true;
bool uniform_sampling = true;
int reduction_strategy = 0;
int primes = 0;

using thread_pool = std::unique_ptr<BS::thread_pool<BS::none>>;

std::vector<int32_t> PRIMES;

// NOTE: PRIMES is still empty during static initialisation of this translation
// unit, and BS::thread_pool treats a count of 0 as "use hardware_concurrency()".
// Sizing the pool from PRIMES here therefore silently gave it as many threads as
// the machine has cores, which made the sampling in weighted_sampling() (and thus
// the whole computation) machine-dependent. The real pool is created in verify()
// once PRIMES is filled -- do not size it from PRIMES here.
thread_pool pool = std::make_unique<BS::thread_pool<BS::none>>(1);
}

/*------------------------------------------------------------------------*/
// All state of a run lives in globals, hence runs are serialized.
static std::mutex run_mutex;

/*------------------------------------------------------------------------*/
static void
apply_options(const Options& opts) {
  verbose = opts.verbose;
  mult_spec = (opts.spec == SpecMode::Mult);
  miter_spec = (opts.spec == SpecMode::Miter);
  talisman::primes = opts.primes;
  talisman::reduction_strategy = opts.reduction;
  talisman::custom_linalg = opts.custom_linalg;
  talisman::uniform_sampling = opts.uniform_sampling;
}

/*------------------------------------------------------------------------*/
// Frees everything a run allocated and restores the initial global state.
static void
reset_run_state() {
  talisman::pool = std::make_unique<BS::thread_pool<BS::none>>(1);
  linpoly_cache_clear();
  clear_detection_cache();
  reset_gates();
  reset_aig_parsing();
  talisman::PRIMES.clear();
  moduli.clear();
  two_pow_NN_imgs.clear();
  reset_randomness();
  reset_statistics();
}

/*------------------------------------------------------------------------*/
// The actual verification, returns true iff guess and prove succeeded.
static bool
verify(const Options& opts) {
  parse_aig(opts.input_file.c_str());

  init_gates(mult_spec ? NN * NN / 4 : 0);
  print_hline();

  linpoly_cache_init();

  // The product of the 8 structured primes is only about 2^128.0009, so for
  // NN = 128 the residues of 2^NN are those of a coefficient of about
  // -2^117.4, which would then be removed as a multiple of 2^NN.
  // Hence NN = 128 needs the 16 primes.
  if(primes == 0 and NN < 128)
    talisman::PRIMES
      = { 65519, 65521, 65537, 65539, 65543, 65551, 65557, 65563 };
  else if(primes == 1 and NN < 128)
    talisman::PRIMES
      = { 91373, 79601, 91583, 127837, 46639, 57719, 111497, 43103 };
  if(primes == 0 and NN >= 128)
    talisman::PRIMES
      = { 65497, 65519, 65521, 65537, 65539, 65543, 65551, 65557,
          65563, 65579, 65581, 65587, 65599, 65609, 65617, 65629 };
  else if(primes == 1 and NN >= 128)
    talisman::PRIMES
      = { 76907,  120539, 91513, 45841,  46399,  96601, 72253,  42743,
          122147, 79231,  73961, 119027, 128659, 49727, 119297, 78989 };

  num_primes = talisman::PRIMES.size();

  // Create the real worker pool now that PRIMES is populated (see the note at
  // the definition of `pool` above).
  talisman::pool
    = std::make_unique<BS::thread_pool<BS::none>>(talisman::PRIMES.size());

  CoeffVec::set_moduli();
  init_two_pow_NN_imgs();

  msg("Primes: %s", talisman::primes == 0 ? "structured" : "random");
  msg("Worker threads: %zu", talisman::pool->get_thread_count());
  msg("Reduction strategy: %s",
      talisman::reduction_strategy == 0   ? "hybrid"
      : talisman::reduction_strategy == 1 ? "linear"
                                          : "nonlinear");
  msg("Using built-in structural detection for adders");
  if(talisman::reduction_strategy < 2) {
    msg("Custom linear algebra: %s",
        talisman::custom_linalg ? "enabled" : "disabled");
    msg("Using %s sampling", talisman::uniform_sampling ? "uniform" : "random");
  }

  // Detection of FA/HA/XOR and FSA nodes via built-in structural detection
  auto abc_start = std::chrono::high_resolution_clock::now();
  try {
    if(verbose >= 2) {
      print_aig_stats();
    }

    // Detect HA/FA/XOR first to populate the cache for get_kcover_nodes
    std::vector<HalfAdder> ha_list;
    std::vector<FullAdder> fa_list;
    // NOTE: xor_list is kept for API compatibility but not used since affine
    // cuts cover it
    std::vector<XorRemain> xor_list;
    if(detect_ha_fa_xor(ha_list, fa_list, xor_list)) {
      if(verbose >= 2) {
        // XOR-remaining count is no longer relevant since we use affine cuts
        // instead
        msg("detected FA count: %zu, HA count: %zu",
            fa_list.size(),
            ha_list.size());
      }

      if(!fa_list.empty()) {
        mark_fahaxor(fa_list, 1);
      }
      if(!ha_list.empty()) {
        mark_fahaxor(ha_list, 0);
      }
      // NOTE: xor_remain processing is commented out because affine cuts
      // already cover the linear relations from xor_remaining gates.
      // if(!xor_list.empty()) {
      //   mark_fahaxor(xor_list);
      // }

      if(opts.trust_abc) {
        trust_linpoly_cache();
      } else {
        check_linpoly_cache(opts.sc_depth, prime);
      }
    } else {
      warning("HA/FA/XOR detection failed. Proceeding without it.");
    }

    // Get K-Cover nodes (reverse_fsa_nodes and boundary_nodes)
    // Reuses the adder pairs cached by detect_ha_fa_xor above.
    std::vector<int> reverse_fsa_nodes;
    std::vector<int> boundary_nodes;
    if(get_kcover_nodes(reverse_fsa_nodes, boundary_nodes)) {
      if(verbose >= 2) {
        msg("detected %zu reverse fsa nodes and %zu boundary nodes.",
            reverse_fsa_nodes.size(),
            boundary_nodes.size());
      }
      mark_fsa(reverse_fsa_nodes, boundary_nodes);
    } else {
      warning("K-Cover detection failed. Proceeding without it.");
    }

    clear_detection_cache();

  } catch(const std::exception& e) {
    warning("Candidate detection error: %s. Proceeding without adder "
            "detection.",
            e.what());
  }
  // testing cache linear polynomial
  if(verbose >= 2) {
    print_linpoly_cache();
  }
  auto abc_end = std::chrono::high_resolution_clock::now();
  abc_time = std::chrono::duration<double>(abc_end - abc_start).count();

  LinPolynomial spec_poly = gen_mult_spec_poly();

  bool result = guess_and_prove(spec_poly, opts.sc_depth);

  if(!result) {
    msg("Guess and Prove failed.");
  } else {
    msg("Guess and Prove succeeded.");
  }

  print_statistics();
  return result;
}

/*------------------------------------------------------------------------*/
namespace talisman {

const char*
version() {
  return VERSION;
}

/*------------------------------------------------------------------------*/
Result
run(const Options& opts, const Hooks& hooks) {
  std::lock_guard<std::mutex> lock(run_mutex);

  reset_run_state();
  apply_options(opts);
  set_sat_solver_factory(hooks.sat_solver_factory);
  set_terminate_callback(hooks.terminate);

  std::optional<SignalGuard> sg;
  if(hooks.install_signal_handlers)
    sg.emplace();

  Result res;
  try {
    res.status = verify(opts) ? Status::Verified : Status::NotVerified;
  } catch(const Interrupted&) {
    msg("interrupted");
    res.status = Status::Interrupted;
  } catch(const Error& e) {
    res.status = Status::Error;
    res.error_code = e.code;
    res.error_message = e.message;
  } catch(const std::exception& e) {
    res.status = Status::Error;
    res.error_code = 255;
    res.error_message = e.what();
  }

  sg.reset();
  set_terminate_callback({});
  set_sat_solver_factory({});
  reset_run_state();
  return res;
}

}// namespace talisman
