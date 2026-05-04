/*------------------------------------------------------------------------*/
/*! \file guess_and_prove.cpp
    \brief contains the functions to prove the guess and prove the guess

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "guess_and_prove.hpp"

#include <cadical.hpp>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <map>
#include <mutex>
#include <thread>
#include <sstream>
#include <limits>

#include "linpoly_cache.hpp"
#include "linpolynomial.hpp"
#include "polynomial.hpp"
#include "subcircuit.hpp"
#include "thread_pool.hpp"

using namespace talisman;

namespace talisman {

using thread_pool = std::unique_ptr<BS::thread_pool<BS::none>>;
extern thread_pool pool;
extern int reduction_strategy;

extern int fsa_size;

extern Gate*
gate(int lit);

// Thread-local accumulators for timing
thread_local double local_prove_time = 0.0;
thread_local double local_repair_time = 0.0;
thread_local double local_guess_time = 0.0;

// Per-thread timing data structure
struct ThreadTimingData {
  double prove_time = 0.0;
  double repair_time = 0.0;
  double guess_time = 0.0;
};

// Map to store per-prime timing data (keyed by prime number)
std::map<int32_t, ThreadTimingData> thread_timing_map;
std::mutex thread_timing_mutex;

/*-------------------------------------------------------------------------*/
// Aggregate thread-local times to global atomics and per-prime map
void aggregate_thread_times(int32_t prime) {
  // Use fetch_add for lock-free atomic addition
  prove_time.fetch_add(local_prove_time, std::memory_order_relaxed);
  repair_time.fetch_add(local_repair_time, std::memory_order_relaxed);
  guess_time.fetch_add(local_guess_time, std::memory_order_relaxed);
  
  // Store in per-prime map
  {
    std::lock_guard<std::mutex> lock(thread_timing_mutex);
    thread_timing_map[prime].prove_time += local_prove_time;
    thread_timing_map[prime].repair_time += local_repair_time;
    thread_timing_map[prime].guess_time += local_guess_time;
  }
  
  // Reset thread-local accumulators
  local_prove_time = 0.0;
  local_repair_time = 0.0;
  local_guess_time = 0.0;
}

/*-------------------------------------------------------------------------*/
// Print individual thread timing statistics
void print_individual_thread_times() {
  std::lock_guard<std::mutex> lock(thread_timing_mutex);
  
  if(thread_timing_map.empty()) {
    msg("No per-thread timing data available");
    return;
  }
  
  // Initialize min/max values
  double min_guess = std::numeric_limits<double>::max();
  double max_guess = 0.0;
  double min_prove = std::numeric_limits<double>::max();
  double max_prove = 0.0;
  double min_repair = std::numeric_limits<double>::max();
  double max_repair = 0.0;
  double min_total = std::numeric_limits<double>::max();
  double max_total = 0.0;
  
  msg("--- Individual Thread Times ---");
  int thread_num = 0;
  for(const auto& [thread_id, timing] : thread_timing_map) {
    double total_time = timing.guess_time + timing.prove_time + timing.repair_time;
    
    std::ostringstream oss;
    oss << thread_id;
    msg("  Thread %d (Prime: %s):", thread_num++, oss.str().c_str());
    msg("    Guess Time:  %.6f seconds", timing.guess_time);
    msg("    Prove Time:  %.6f seconds", timing.prove_time);
    msg("    Repair Time: %.6f seconds", timing.repair_time);
    msg("    Total Time:  %.6f seconds", total_time);
    
    // Update min/max
    min_guess = std::min(min_guess, timing.guess_time);
    max_guess = std::max(max_guess, timing.guess_time);
    min_prove = std::min(min_prove, timing.prove_time);
    max_prove = std::max(max_prove, timing.prove_time);
    min_repair = std::min(min_repair, timing.repair_time);
    max_repair = std::max(max_repair, timing.repair_time);
    min_total = std::min(min_total, total_time);
    max_total = std::max(max_total, total_time);
  }
  
  msg("");
  msg("--- Min/Max Thread Times ---");
  msg("  Guess Time:  min=%.6f, max=%.6f seconds", min_guess, max_guess);
  msg("  Prove Time:  min=%.6f, max=%.6f seconds", min_prove, max_prove);
  msg("  Repair Time: min=%.6f, max=%.6f seconds", min_repair, max_repair);
  msg("");
}

/*-------------------------------------------------------------------------*/

LinPolynomial
gen_output_poly(Gate* g) {
  LinPolynomial lp;
  if(!g)
    return lp;
  // only one child (second is true)
  if(g->get_ch1() == 1) {
    lp.add_term(g, 1);
    Gate* c0 = gate(g->get_ch0());

    if(g->get_ch0() & 1) {
      lp.add_term(c0, 1);
      lp.add_term(nullptr, -1);
    } else {
      lp.add_term(c0, -1);
    }
    return lp;
  }
}

/*-------------------------------------------------------------------------*/
bool
try_eliminate_linear_output(Gate* g, LinPolynomial& tmp) {
  if(g->get_ch1() != 1)
    return false;

  if(verbose > 0) {
    msg("having single child node %s", g->get_var_name().c_str());
  }

  LinPolynomial out = gen_output_poly(g);
  if(out.get_size() == 0)
    return false;

  if(verbose > 0) {
    msg("Subtracting output polynomial:");
    out.print();
  }

  
  tmp.reduce_by(out);


  return true;
}
/*-------------------------------------------------------------------------*/
thread_local std::vector<LinPolynomial> guesses;
thread_local std::vector<std::vector<int>> repair_samples;
thread_local std::unordered_set<Gate*> already_proven;
bool
prove_and_repair_guesses(Subcircuit& sc,
                         nmod_mat_t mat,
                         std::vector<LinPolynomial>& input) {
  bool wrong = true;
  bool candidate = true;
  int num_repair = 0;

  guesses = std::move(input);
  input.clear();
  repair_samples.clear();
  already_proven.clear();

  int prime = mat->mod.n;

  Gate* root = nullptr;
  if(guesses.size() > 0)
    root = guesses[0].get_leading_term().first;

  while(wrong) {
    repair_samples.clear();

    auto prove_start = std::chrono::high_resolution_clock::now();

    CaDiCaL::Solver* solver = new CaDiCaL::Solver;
    int32_t fresh_var = encode_subcircuit(sc, solver);
    wrong = false;

    for(auto& g : guesses) {
      Gate* lm = g.get_leading_term().first;
      // don't prove guess if we already have correct poly
      if(already_proven.contains(lm) or lm->get_gatepoly())
        continue;

      int coeff_sum = g.minimize_coeffs(prime);
      if(lm != root
         && (coeff_sum > 10000 || g.get_leading_term().second > 1000)) {
        if(verbose > 2)
          msg("Skipping prove for gate %s: coeffs too large",
              lm->get_var_name().c_str());
        continue;
      }

      num_proves++;

      bool guess_is_wrong
        = prove_guess(sc, g, solver, repair_samples, prime, fresh_var);
      wrong = wrong or guess_is_wrong;

      if(not guess_is_wrong) {
        num_correct_guesses++;
        already_proven.insert(lm);
        input.push_back(g);
      }
      if(repair_samples.size() >= 200)
        break;
    }
    delete solver;
    auto prove_end = std::chrono::high_resolution_clock::now();
    local_prove_time
      += std::chrono::duration<double>(prove_end - prove_start).count();

    auto repair_start = std::chrono::high_resolution_clock::now();
    if(wrong) {
      num_repair_loop++;
      num_repair++;
      std::tie(candidate, guesses) = repair_guess(sc, repair_samples, mat);
    }
    auto repair_end = std::chrono::high_resolution_clock::now();
    local_repair_time
      += std::chrono::duration<double>(repair_end - repair_start).count();

    if(!candidate) {
      msg("No candidates found after repair");
      input.clear();
      return false;
    }
  }

  if(num_repair > max_num_repair_loop)
    max_num_repair_loop = num_repair;

  return true;
}
/*-------------------------------------------------------------------------*/
std::unordered_map<Gate*, std::vector<LinPolynomial>> poly_map;
std::vector<LinPolynomial> out;

std::vector<LinPolynomial>
combine_multimod_polies(
  std::vector<std::vector<LinPolynomial>>& modular_guesses) {

  poly_map.clear();
  for(auto& guesses_p : modular_guesses)
    for(auto g : guesses_p)
      poly_map[g.get_leading_term().first].push_back(g);

  out.clear();
  out.reserve(poly_map.size());

  for(auto& [lm, polies] : poly_map) {
    LinPolynomial p;
    for(int i = 0; i < (int)polies.size(); i++) {
      auto poly = polies[i];
      for(auto& term : poly.get_terms())
        p.add_term(term.first, term.second, i);
    }
    if(p.multimod_terms.begin()->second.is_nonzero())
      out.push_back(p);
  }

  return out;
}
/*-------------------------------------------------------------------------*/
std::unordered_set<Gate*> reachable_nodes;
std::queue<Gate*> to_visit;
std::vector<Gate*> filtered_nodes;
std::vector<Gate*> filtered_inputs;
std::unordered_set<Gate*> seen;

static void
filter_large_subcircuit(Gate* leading, Subcircuit& sc) {
  // Filter subcircuit by tracing paths from leading node to inputs
  std::unordered_set<Gate*> original_inputs(sc.inputs.begin(), sc.inputs.end());
  reachable_nodes.clear();

  // Always start from leading node and trace backward through data dependencies
  to_visit.push(leading);
  reachable_nodes.insert(leading);

  // Also include an AIG output with cone_id exactly one less than the leading
  // node
  int root_cone = leading->get_cone_id();
  if(verbose > 0) {
    msg("Leading node %s has cone_id=%d",
        leading->get_var_name().c_str(),
        root_cone);
  }

  for(Gate* g : sc.nodes) {
    if(g->get_cone_id() >= root_cone - 2 && g->get_cone_id() <= root_cone +1
       && g->get_level() <= leading->get_level()) {
      if(reachable_nodes.find(g) == reachable_nodes.end()) {
        to_visit.push(g);
        reachable_nodes.insert(g);
      }
    }
  }

  // BFS backward through circuit from leading node to inputs
  while(!to_visit.empty()) {
    Gate* current = to_visit.front();
    to_visit.pop();

    // Get children (predecessor nodes in data flow)
    unsigned ch0_lit = current->get_ch0();
    unsigned ch1_lit = current->get_ch1();

    // Process ch0 if it's not a constant
    if(ch0_lit > 1) {
      Gate* child0 = gate(ch0_lit);
      if(child0 && !original_inputs.count(child0)
         && reachable_nodes.find(child0) == reachable_nodes.end()) {
        reachable_nodes.insert(child0);
        // Continue traversing unless it's an input
        if(!child0->get_input()) {
          to_visit.push(child0);
        }
      }
    }

    // Process ch1 if it's not a constant and not just a literal
    if(ch1_lit > 1) {
      Gate* child1 = gate(ch1_lit);
      if(child1 && !original_inputs.count(child1)
         && reachable_nodes.find(child1) == reachable_nodes.end()) {
        reachable_nodes.insert(child1);
        // Continue traversing unless it's an input
        if(!child1->get_input()) {
          to_visit.push(child1);
        }
      }
    }
  }

  // Filter nodes to only include those reachable from leading node to inputs
  filtered_nodes.clear();
  filtered_inputs.clear();
  for(Gate* g : sc.nodes) {
    if(reachable_nodes.find(g) != reachable_nodes.end()) {
      filtered_nodes.push_back(g);
      if(g->get_ch0() > 1) {
        Gate* c0 = gate(g->get_ch0());
        if(reachable_nodes.find(c0) == reachable_nodes.end())
          filtered_inputs.push_back(c0);
      }
      if(g->get_ch1() > 1) {
        Gate* c1 = gate(g->get_ch1());
        if(reachable_nodes.find(c1) == reachable_nodes.end())
          filtered_inputs.push_back(c1);
      }
    }
  }

  seen.clear();
  auto it = filtered_inputs.begin();
  for(auto& x : filtered_inputs) {
    if(seen.insert(x).second) {
      *it++ = x;
    }
  }
  filtered_inputs.erase(it, filtered_inputs.end());

  size_t original_size = sc.inputs.size() + sc.nodes.size();
  size_t filtered_size = filtered_inputs.size() + filtered_nodes.size();
  if(filtered_size > 0.6 * original_size)
    return;

  sc.inputs = std::move(filtered_inputs);
  sc.nodes = std::move(filtered_nodes);

  if(verbose > 0) {
    msg("Filtered subcircuit from %zu to %zu gates (keeping nodes on paths "
        "from leading node to inputs)",
        original_size,
        filtered_size);
  }
}

/*-------------------------------------------------------------------------*/

bool
get_guesses_and_prove(Gate* leading,
                      Subcircuit& sc,
                      std::vector<LinPolynomial>& correct_guesses) {

  // Filter subcircuit if FSA size exceeds threshold and leading gate is part of
  // FSA Once started, always filter as long as leading gate is part of FSA
  

  if(fsa_size > 3000 && leading->is_fsa() && get_size(sc) > 500) {
    if(verbose > 0) {
      msg("FSA size %zu exceeds threshold and leading gate is FSA, filtering "
          "subcircuit",
          fsa_size);
    }
    filter_large_subcircuit(leading, sc);
    
  }

  size_t samples = std::min(3 * get_size(sc), 20'000UL);

  // sample circuit
  nmod_mat_t mat;
  auto sample_start = std::chrono::high_resolution_clock::now();
  sample_and_eval_subcircuit(sc, samples, mat);
    auto sample_end = std::chrono::high_resolution_clock::now();
  sample_time
    += std::chrono::duration<double>(sample_end - sample_start).count();

  // Check if sample_and_eval_subcircuit successfully initialized the matrix
  // Early returns in sample_and_eval_subcircuit leave mat uninitialized
  if(!mat || mat->r == 0 || mat->c == 0) {
    return false;
  }

  if(verbose > 0)
    msg("Starting linear algebra for matrix of size %d x %d", mat->r, mat->c);

  std::vector<std::vector<LinPolynomial>> modular_guesses(PRIMES.size());
  for(size_t i = 0; i < PRIMES.size(); i++) {

    auto p = PRIMES[i];

    auto task = [&sc, &mat, &modular_guesses, p, i] {
      auto guess_start = std::chrono::high_resolution_clock::now();
      nmod_mat_t mat_p;
      nmod_mat_init_set(mat_p, mat);
      nmod_mat_set_mod(mat_p, p);

      bool candidate;
      std::tie(candidate, modular_guesses[i]) = guess_naive(sc, mat_p);

      if(verbose > 0) {
        msg("Modular guessing done for prime %d, found %zu candidates",
            p,
            modular_guesses[i].size());
      }
        
      num_guesses += modular_guesses[i].size();
      auto guess_end = std::chrono::high_resolution_clock::now();
      
      local_guess_time
        += std::chrono::duration<double>(guess_end - guess_start).count();

      if(candidate)
        prove_and_repair_guesses(sc, mat_p, modular_guesses[i]);
      else
        modular_guesses[i].clear();

      nmod_mat_clear(mat_p);
      
      // Aggregate thread-local times to global atomics
      aggregate_thread_times(p);
    };
    pool->detach_task(task);
  }
  pool->wait();

  // if any of entries is empty, return false
  if(std::any_of(modular_guesses.begin(), modular_guesses.end(), [](auto& v) {
       return v.empty();
     })) {
    nmod_mat_clear(mat);
    return false;
  }


  auto lin_start = std::chrono::high_resolution_clock::now();
  correct_guesses = combine_multimod_polies(modular_guesses);
  auto lin_end = std::chrono::high_resolution_clock::now();
  lin_time
      += std::chrono::duration<double>(lin_end - lin_start).count();
  
  nmod_mat_clear(mat);

  return true;
}
/*-------------------------------------------------------------------------*/
bool
reduce_by_correct_guesses(LinPolynomial& tmp,
                          std::vector<LinPolynomial>& correct) {

  for(auto& g : correct) {
    auto lm = g.get_leading_term().first;
    if(tmp.get_leading_term().first->get_level() != lm->get_level()) {
      if(!lm->get_gatepoly()) {
        lm->set_gatepoly(g);
      }
      continue;
    }

    if(verbose > 0) {
      msg("Reduce by correct guess:");
      g.print_multimod();
    }

   


    
    tmp.reduce_by(g);

    if(tmp.get_size() == 0)
      return true;
  }
  return false;
}

/*-------------------------------------------------------------------------*/
bool
increase_subcircuit_or_fail(int& sc_depth, Gate* leading, Subcircuit& sc) {
  if(verbose > 0) {
    msg("Trying to increase subcircuit depth from %d for leading gate %s "
        "(level=%d, dist=%d)",
        sc_depth,
        leading->get_var_name().c_str(),
        leading->get_level(),
        leading->get_dist());
  }
  // If requested depth is strictly less than the distance to inputs, we
  // can safely increase the depth.
  if(sc_depth < (int)leading->get_dist()) {
    sc_depth += sc_depth < 4 ? 2 : 1;
    return true;
  }

  // If the requested depth equals the distance, prefer to try a
  // neighbor-COI expansion first (keeps the same depth) when the current
  // subcircuit is smaller than the gate's level. This allows including
  // nodes from one smaller COI as an alternative root while preserving
  // the distance/level budget.
  if(sc_depth >= (int)leading->get_dist()) {
    if(sc.nodes.size() < static_cast<size_t>(leading->get_level()) - NN) {
      if(verbose > 0) {
        msg("increase_subcircuit_or_fail: keeping depth=%d and allowing "
            "neighbor-COI expansion for %s",
            sc_depth,
            leading->get_var_name().c_str());
      }
      return true;
    }
    // otherwise increase depth as a fallback
  }

  return false;
}

/*-------------------------------------------------------------------------*/
bool
try_reduce_by_stored_poly(Gate* lm, LinPolynomial& tmp) {
  if(!lm->get_gatepoly())
    return false;

  LinPolynomial stored_poly = *(lm->get_gatepoly());
  if(verbose > 0) {
    msg("Found stored polynomial for gate %s, trying to reduce.",
        lm->get_var_name().c_str());
    stored_poly.print_multimod();
  }


  tmp.reduce_by(stored_poly);

  return true;
}
/*-------------------------------------------------------------------------*/
void
find_carries() {
  for(unsigned j = num_gates - 1; j > 0; j--) {
    Gate* n = gates[j];
    if(n->get_input())
      continue;
    if(n->is_extension())
      continue;

    n->set_carry_gate(0);
    const auto& parents = n->get_parents();
    for(unsigned parent_lit : parents) {
      Gate* n_parent = gate(parent_lit);
      if(n_parent && n_parent->get_cone_id() > n->get_cone_id()) {
        n->inc_carry_gate();
      }
    }
  }
}

bool
try_nonlinear_rewriting(LinPolynomial& tmp) {

  if(verbose > 0)
    msg("Attempting nonlinear rewriting fallback");

    

  if(tmp.get_leading_term().first->get_dist() > 6) {
    // reassign_levels_by_cone_order();
    auto simplify_start = std::chrono::high_resolution_clock::now();
    find_carries();
    decompose_all_gate_polynomials(tmp.get_leading_term().first);
    auto simplify_end = std::chrono::high_resolution_clock::now();
    simplify_time += std::chrono::duration<double>(simplify_end - simplify_start).count();
  }
    

  // print_cones();

  // adjust_cones_by_merging();

  translate_all_stored_gate_polys();

  if(tmp.get_size() == 0) {
    if(verbose > 0) {
      msg("Polynomial is already zero, no rewriting needed");
    }
    return true;
  }

  // Get the leading gate from the LinPolynomial
  Gate* leading_gate = tmp.get_leading_term().first;
  if(!leading_gate) {
    if(verbose > 0) {
      msg("No leading gate in polynomial for nonlinear rewriting");
    }
    return false;
  }

  // Convert LinPolynomial to Polynomial
  Polynomial poly_tmp = Polynomial::from_linpoly(tmp);

  if(verbose > 0) {
    msg("Original polynomial:");
    poly_tmp.print_multimod();
  }

  // Keep reducing by gate polynomials until leading gate is an input or
  // polynomial is zero
  while(!poly_tmp.is_zero()) {

    if(poly_tmp.check_leading_coeff())
      continue;

    // Get the leading monomial
    const auto& leading_term = poly_tmp.get_lm();

    // Find the leading gate (highest level variable in leading monomial)
    Gate* leading_var = nullptr;
    for(Gate* g : leading_term.vars) {
      if(g) {
        leading_var = g;
        break;// Already sorted by level, so first is highest
      }
    }

    if(!leading_var) {
      // Reached a constant, polynomial is reduced
      break;
    }

    // Check if leading variable is an input
    if(leading_var->get_input()) {
      if(verbose > 0) {
        msg("Reached input variable %s, nonlinear reduction complete",
            leading_var->get_var_name().c_str());
      }
      break;
    }

    if(verbose > 0) {
      msg("Reducing by gate polynomial for %s (level %d, cone %d, dist %d)",
          leading_var->get_var_name().c_str(),
          leading_var->get_level(),
          leading_var->get_cone_id(),
          leading_var->get_dist());
    }

    Polynomial gate_poly;

    if(leading_var->has_gatepoly() && leading_var->get_faha_type() == 2) {
      if(verbose > 0) {
        msg("Using stored gate polynomial for %s",
            leading_var->get_var_name().c_str());
      }
      // Use stored gate polynomial if available
      Polynomial translated
        = Polynomial::from_linpoly(*(leading_var->get_gatepoly()));

      gate_poly = translated;
    }

    if(gate_poly.is_zero() && leading_var->has_gatepoly_nonlinear()) {
      if(verbose > 0) {
        msg("Using stored nonlinear gate polynomial for %s",
            leading_var->get_var_name().c_str());
      }
      // Use stored nonlinear gate polynomial if available
      gate_poly = *(leading_var->get_gatepoly_nonlinear());
    }

    if(gate_poly.is_zero()) {
      gate_poly = Polynomial::get_gate_poly(leading_var);
    }

    if(gate_poly.get_size() == 1) {
      gate_poly.print_multimod();
      throw std::runtime_error("Gate polynomial has size 1, cannot reduce");
    }

    if(verbose > 0) {
      msg("Gate polynomial (multimod):");
      gate_poly.print_multimod();
    }

    // Reduce the polynomial by the gate polynomial
    poly_tmp.reduce_by(gate_poly);
    // msg("size of poly after reduction: %u", poly_tmp.get_size());

    if(verbose >= 2) {
      msg("After reduction:");
      poly_tmp.print_multimod();
      msg("");
    }
  }

  msg("Nonlinear reduction complete:");
  poly_tmp.print_multimod();

  // msg("Size after nonlinear reduction: %zu", poly_tmp.get_size());

  // Successfully reduced using nonlinear approach
  // We keep the result as a polynomial (not converting back to linear)
  // The caller will need to handle this appropriately
  return poly_tmp.is_zero();
}
/*-------------------------------------------------------------------------*/

std::vector<LinPolynomial> correct_guesses;
bool
compute_linear_polynomials(Subcircuit& sc, Gate* leading) {

  correct_guesses.clear();

  bool res = false;

  // small subcircuit -> compute full truth table
  if(sc.inputs.size() <= 12){

    auto guess_start = std::chrono::high_resolution_clock::now();
    res = linear_polies_from_truth_table(sc, correct_guesses);
    auto guess_end = std::chrono::high_resolution_clock::now();
    guess_noparallel_time
      += (std::chrono::duration<double>(guess_end - guess_start).count());
  }
    
  else
    // large subcircuit -> do guess & prove
    res = get_guesses_and_prove(leading, sc, correct_guesses);

  if(!res)
    return false;

  // Check that at least one correct guess has leading_gate as its leading
  // variable
  bool found_matching_poly = false;
  for(const auto& guess : correct_guesses) {
    auto [guess_leading, guess_coeff] = guess.get_leading_term();
    if(guess_leading == leading) {
      found_matching_poly = true;
      break;
    }
  }

  if(!found_matching_poly) {
    if(verbose > 0) {
      msg("No correct guess found with leading variable %s",
          leading->get_var_name().c_str());
    }
    return false;
  }
  return true;
}

/*-------------------------------------------------------------------------*/
bool
guess_and_prove(LinPolynomial& poly, int sc_depth_init) {
  LinPolynomial tmp = poly;
  int sc_depth = sc_depth_init;
  Gate* leading = tmp.get_leading_term().first;
  bool after_fsa = false;

  Gate* tmp_leading = nullptr;
  bool new_leading = false;

  if(reduction_strategy >= 2) {
    auto nonlin_start = std::chrono::high_resolution_clock::now();
    tmp.multimod = true;
    bool nonlin_res = try_nonlinear_rewriting(tmp);
    auto nonlin_end = std::chrono::high_resolution_clock::now();
    nonlin_time
      += std::chrono::duration<double>(nonlin_end - nonlin_start).count();
    return nonlin_res;
  }

  while(tmp.get_size() > 0) {
    tmp.multimod = true;

    // Check if leading coeff is multiple 
    auto lin_start = std::chrono::high_resolution_clock::now();
    bool leading_coeff_large = tmp.check_leading_coeff();
    auto lin_end = std::chrono::high_resolution_clock::now();
    lin_time
      += std::chrono::duration<double>(lin_end - lin_start).count();
    if(leading_coeff_large) continue;

  
    
    leading = tmp.get_leading_term().first;
    
    // test if leading monmial is constant
    if(leading == nullptr) {
      msg("Polynomial reduced to constant only.");
      tmp.print_multimod();
      return false;
    }

    // Reset sc depth if new leading term was found 
    if(leading != tmp_leading) {
      new_leading = true;
      tmp_leading = leading;
      sc_depth = sc_depth_init;
    } else {
      new_leading = false;
    }

     // Phase 1: trivial linear output elimination
    lin_start = std::chrono::high_resolution_clock::now();
    bool elim_lin_output = try_eliminate_linear_output(leading, tmp);
    lin_end = std::chrono::high_resolution_clock::now();
    lin_time
      += std::chrono::duration<double>(lin_end - lin_start).count();
    if(elim_lin_output) continue; 

  
    // Phase 2: check if a correct polynomial is stored for this gate
    lin_start = std::chrono::high_resolution_clock::now();
    bool elim_stored_output = try_reduce_by_stored_poly(leading, tmp);
    lin_end = std::chrono::high_resolution_clock::now();
    lin_time
      += std::chrono::duration<double>(lin_end - lin_start).count();
    if(elim_stored_output) continue; 

    
    if(fsa_size > 10000 && !leading->is_fsa() && !after_fsa) {
      reset_levels_to_orig();
      after_fsa = true;
    }

    // Phase 3: build subcircuit
    auto sc_start = std::chrono::high_resolution_clock::now();
    Subcircuit sc = get_subcircuit(leading, sc_depth, new_leading);
    auto sc_end = std::chrono::high_resolution_clock::now();
    extract_sc_time
          += std::chrono::duration<double>(sc_end - sc_start).count();
   


    bool linear_rewriting = compute_linear_polynomials(sc, leading);

    // Phase 4: get linear polynomials
    if(!linear_rewriting) {
      if(sc_depth < 5 or reduction_strategy == 1) {
        if(increase_subcircuit_or_fail(sc_depth, leading, sc))
          continue;// Successfully increased, try again with larger subcircuit
      }

      // Cannot increase subcircuit -> use nonlinear fallback
      if(tmp.get_size() > 0) {
        if(verbose > 0) {
          msg("Fallback, trying nonlinear rewriting fallback after reduction "
              "failure",
              leading->get_var_name().c_str(),
              leading->get_dist());
        }
        auto nonlin_start = std::chrono::high_resolution_clock::now();
        bool nonlin_res = try_nonlinear_rewriting(tmp);
        auto nonlin_end = std::chrono::high_resolution_clock::now();
        nonlin_time
          += std::chrono::duration<double>(nonlin_end - nonlin_start).count();
        return nonlin_res;
      }
    }

    // Phase 5: reduce polynomial
   
    lin_start = std::chrono::high_resolution_clock::now();
    bool reduced_guesses = reduce_by_correct_guesses(tmp, correct_guesses);
    lin_end = std::chrono::high_resolution_clock::now();
    lin_time
      += std::chrono::duration<double>(lin_end - lin_start).count();
    if(reduced_guesses) return true; 

  }

  return tmp.get_size() == 0;
}
}// namespace talisman
