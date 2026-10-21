/*------------------------------------------------------------------------*/
/*! \file prove.cpp
    \brief Contains functions necessary to prove the guess.

    Part of TalisMan
    Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "prove.hpp"
#include "fast_mod.hpp"

#include <atomic>
#include <bit>
#include <cadical.hpp>
#include <unordered_map>

#include "linpoly_cache.hpp"
#include "pblib/pb2cnf.h"
#include "thread_pool.hpp"
#include <iostream>

using namespace talisman;

namespace talisman {

thread_local std::vector<int> sat_model;
static std::vector<int>
get_model(CaDiCaL::Solver* solver, Subcircuit& sc) {
  sat_model.clear();
  sat_model.reserve(get_size(sc) + 1);

  // constant column
  sat_model.push_back(1);

  for(auto it = sc.inputs.rbegin(); it != sc.inputs.rend(); it++) {
    int var = (*it)->get_var() >> 1;
    int val = solver->val(var) > 0;
    sat_model.push_back(val);
  }

  for(auto it = sc.nodes.rbegin(); it != sc.nodes.rend(); it++) {
    int var = (*it)->get_var() >> 1;
    int val = solver->val(var) > 0;
    sat_model.push_back(val);
  }

  return sat_model;
}

//------------------------------------------------------------------------
static int64_t
evaluate(CaDiCaL::Solver* solver, LinPolynomial& poly) {
  int64_t result = 0;
  for(const auto& term : poly.get_terms()) {
    Gate* gate = term.first;
    int32_t coeff = term.second;

    if(gate == nullptr) {
      if(__builtin_add_overflow(result, coeff, &result))
        die(7, "addition overflow in evaluate");

    } else {
      int var = gate->get_var() >> 1;
      bool val = solver->val(var) > 0;
      if(val)
        if(__builtin_add_overflow(result, coeff, &result))
          die(7, "addition overflow in evaluate");
    }
  }

  return result;
}
//------------------------------------------------------------------------
thread_local std::vector<std::vector<int32_t>> formula_geq;
thread_local std::vector<std::vector<int32_t>> formula_leq;

static void
poly_encoding(std::vector<std::vector<int32_t>>& poly_clauses,
              const std::vector<int64_t>& weights,
              const std::vector<int32_t>& literals,
              int32_t constant,
              int prime,
              int32_t& fresh_var,
              int factor) {
  // Use a single global PB2CNF instance protected by mutex.
  // Avoid thread_local to prevent use-after-free during TLS cleanup
  // when complex objects with shared_ptr members are destroyed.
  static PB2CNF pb2cnf;
  static std::mutex pb2cnf_mutex;

  poly_clauses.clear();
  formula_geq.clear();
  formula_leq.clear();

  // Constraint: polynomial > 0
  {
    std::lock_guard<std::mutex> lock(pb2cnf_mutex);

    fresh_var = pb2cnf.encodeGeq(weights,
                                 literals,
                                 factor * prime - constant + 1,
                                 formula_geq,
                                 fresh_var)
                + 1;

    // Constraint: polynomial < 0
    fresh_var = pb2cnf.encodeLeq(weights,
                                 literals,
                                 factor * prime - constant - 1,
                                 formula_leq,
                                 fresh_var)
                + 1;
  }

  //------------------------------------------------------------------------
  // Combine constraints: (a_k ∨ (poly > k*prime)) ∧ (¬a_k ∨ (poly < k*prime))
  //------------------------------------------------------------------------
  int32_t assumption_var = fresh_var++;

  for(auto& clause : formula_geq)
    clause.push_back(assumption_var);

  for(auto& clause : formula_leq)
    clause.push_back(-assumption_var);

  poly_clauses.insert(poly_clauses.end(),
                      std::make_move_iterator(formula_geq.begin()),
                      std::make_move_iterator(formula_geq.end()));
  poly_clauses.insert(poly_clauses.end(),
                      std::make_move_iterator(formula_leq.begin()),
                      std::make_move_iterator(formula_leq.end()));
}
//------------------------------------------------------------------------
int32_t
encode_subcircuit(Subcircuit& sc, CaDiCaL::Solver* solver) {
  //------------------------------------------------------------------------
  // Encode subcircuit CNF into solver
  //------------------------------------------------------------------------
  int32_t fresh_var = num_gates + 1;

  for(const auto& gate : sc.nodes) {
    // if (gate->get_level() > poly.get_leading_term().first->get_level())
    // {
    //   continue; // Skip gates that are not relevant for the polynomial
    // }
    const std::vector<int> cnf = gate->get_cnf_encoding();
    for(int clause : cnf) {
      solver->add(clause);
    }
  }

  return fresh_var;
}
//------------------------------------------------------------------------
thread_local std::vector<std::vector<int32_t>> poly_clauses;
thread_local std::vector<int64_t> weights;
thread_local std::vector<int32_t> literals;

bool
prove_guess_SAT(Subcircuit& sc,
                LinPolynomial& poly,
                CaDiCaL::Solver* solver,
                std::vector<std::vector<int>>& repair_samples,
                int prime,
                int32_t& fresh_var,
                std::vector<int> factors) {

  //------------------------------------------------------------------------
  // Encode linear polynomial in CNF
  //------------------------------------------------------------------------

  weights = poly.get_coeffs_64();
  literals = poly.get_literals();
  const int32_t constant = poly.get_constant();

  for(int factor : factors)
    poly_encoding(
      poly_clauses, weights, literals, constant, prime, fresh_var, factor);

  int selector_var = fresh_var++;
  for(const auto& clause : poly_clauses) {
    for(int lit : clause)
      solver->add(lit);
    solver->add(selector_var);
    solver->add(0);
  }

  //------------------------------------------------------------------------
  // Solve
  //------------------------------------------------------------------------
  solver->assume(-selector_var);
  int result = solver->solve();

  bool res = false;
  while(result == 10) {
    int64_t eval = evaluate(solver, poly);
    // poly.print();
    // std::cout << "SATISFIABLE. Evaluation: " << eval << std::endl;

    if(mod(eval, prime, 0) != 0) {
      if(verbose > 0) {
        std::cout
          << "Polynomial does not evaluate to zero modulo, repair step needed "
          << prime << std::endl;
      }
      repair_samples.push_back(get_model(solver, sc));
      res = true;
      break;
    } else {
      // std::cout << "Polynomial evaluates to zero modulo " << prime <<
      // std::endl; divide eval/prime to get factor
      int64_t factor = eval / prime;
      poly_encoding(
        poly_clauses, weights, literals, constant, prime, fresh_var, factor);
      for(const auto& clause : poly_clauses) {
        for(int lit : clause)
          solver->add(lit);
        solver->add(selector_var);
        solver->add(0);
      }
      // time the subsequent SAT solve
      solver->assume(-selector_var);
      result = solver->solve();
    }
  }

  // poly.print();
  solver->add(selector_var);
  solver->add(0);
  return res;
}

//------------------------------------------------------------------------
bool
prove_guess(Subcircuit& sc,
            LinPolynomial& guess,
            CaDiCaL::Solver* solver,
            std::vector<std::vector<int>>& repair_samples,
            int prime,
            int32_t& fresh_var,
            std::vector<int> factors) {

  if(verbose > 0) {
    std::cout << "proving guess: ";
    guess.print();
  }
  return prove_guess_SAT(
    sc, guess, solver, repair_samples, prime, fresh_var, factors);
}

/*-------------------------------------------------------------------------*/

// Structure to hold task data for parallel linpoly checking
struct LinpolyCheckTask {
  const Gate* gate;
  LinPolynomial poly;
  size_t poly_idx;
};

// Structure to hold result of a single check
struct LinpolyCheckResult {
  const Gate* gate;
  LinPolynomial poly;
  size_t poly_idx;
  bool success;
};

bool
check_linpoly_cache(int sc_depth, int prime, std::vector<int> factors) {
  // Collect all tasks first
  std::vector<LinpolyCheckTask> tasks;

  linpoly_cache_for_each(
    [&](const Gate* g, const std::vector<LinPolynomial>& cached) {
      if(!g)
        return;

      for(size_t j = 0; j < cached.size(); ++j) {
        LinPolynomial lp = cached[j];
        Gate* lead = lp.get_leading_term().first;

        if(!lead) {
          if(verbose >= 1)
            std::cout
              << "[linpoly_check] cached poly has no leading gate, skipping\n";
          continue;
        }

        // Skip XOR-substitute cases for now (they need special handling)
        if(lead->get_xor_root_id() != -1) {
          continue;
        }

        tasks.push_back({ g, lp, j });
      }
    });

  if(tasks.empty()) {
    if(verbose >= 1)
      std::cout << "[linpoly_check] No polynomials to check.\n";
    return true;
  }

  if(verbose >= 1)
    std::cout << "[linpoly_check] Checking " << tasks.size()
              << " polynomials in parallel...\n";

  // Results vector (one per task)
  std::vector<LinpolyCheckResult> results(tasks.size());

  // Atomic counters for progress
  std::atomic<int> success_count{ 0 };
  std::atomic<int> failure_count{ 0 };

  // Create a local thread pool for this operation. Sized like the global pool
  // (one worker per prime) rather than from hardware_concurrency(), which would
  // stack a second core-count-sized pool on top of the global one and make the
  // amount of parallelism machine-dependent.
  BS::thread_pool<BS::none> local_pool(PRIMES.size());

  // Submit all tasks to the thread pool
  for(size_t i = 0; i < tasks.size(); ++i) {
    auto task = [i,
                 &tasks,
                 &results,
                 &success_count,
                 &failure_count,
                 sc_depth,
                 prime,
                 factors]() {
      const LinpolyCheckTask& t = tasks[i];
      LinpolyCheckResult& r = results[i];

      r.gate = t.gate;
      r.poly = t.poly;
      r.poly_idx = t.poly_idx;

      Gate* lead = t.poly.get_leading_term().first;

      // Build subcircuit
      Subcircuit sc = get_subcircuit(lead, sc_depth, 1);

      // Create solver and encode
      CaDiCaL::Solver* solver = new CaDiCaL::Solver;
      int32_t fresh_var = encode_subcircuit(sc, solver);

      // Local repair samples (not used but required by API)
      std::vector<std::vector<int>> local_repair;

      // Make a mutable copy of the polynomial
      LinPolynomial lp_copy = t.poly;
      abc_polies++;
      bool wrong = prove_guess(
        sc, lp_copy, solver, local_repair, prime, fresh_var, factors);
      delete solver;

      r.success = !wrong;

      if(wrong) {
        failure_count.fetch_add(1, std::memory_order_relaxed);
      } else {
        success_count.fetch_add(1, std::memory_order_relaxed);
      }
    };

    local_pool.detach_task(task);
  }

  // Wait for all tasks to complete
  local_pool.wait();

  // Process results: attach successful polynomials to gates
  // This must be done serially to avoid race conditions on gate modification
  for(const auto& r : results) {
    if(r.success) {
      if(verbose >= 1)
        std::cout << "[linpoly_check] OK for gate " << r.gate->get_var_name()
                  << " poly[" << r.poly_idx << "]" << std::endl;

      // Attach verified polynomial to the gate
      Gate* ng = const_cast<Gate*>(r.gate);
      ng->set_gatepoly(r.poly);
    } else {
      Gate* lead = r.poly.get_leading_term().first;
      if(verbose >= 1)
        std::cout << "[linpoly_check] FAILED for gate "
                  << r.gate->get_var_name() << " poly[" << r.poly_idx << "]"
                  << " (lead=" << (lead ? lead->get_var_name() : "null") << ")"
                  << std::endl;
      warning("Linpoly cache check failed for gate %s poly[%zu]",
              r.gate->get_var_name().c_str(),
              r.poly_idx);
    }
  }

  int final_success = success_count.load();
  int final_failure = failure_count.load();

  if(verbose >= 1)
    std::cout << "[linpoly_check] Summary: " << final_success << " succeeded, "
              << final_failure << " failed." << std::endl;

  return final_failure == 0;
}

//------------------------------------------------------------------------
bool
trust_linpoly_cache() {
  // Iterate over the cache and directly attach cached polynomials to their
  // gates without performing SAT/MIP proving. This trusts earlier ABC
  // processing.
  linpoly_cache_for_each(
    [&](const Gate* g, const std::vector<LinPolynomial>& cached) {
      if(!g)
        return;

      if(cached.empty())
        return;

      Gate* ng = const_cast<Gate*>(g);
      // Attach the first cached polynomial by default, and if multiple
      // polynomials exist, attach them in order by calling set_gatepoly
      // repeatedly (assumes set_gatepoly stores the polynomial
      // appropriately).
      for(const auto& lp : cached) {
        abc_polies++;
        ng->set_gatepoly(lp);
      }
    });

  return true;
}

}// namespace talisman
