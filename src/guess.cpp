/*------------------------------------------------------------------------*/
/*! \file guess.cpp
    \brief Contains functions necessary to guess.

    Part of TalisMan
    Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "guess.hpp"
#include "linear_algebra.hpp"
#include "nmod.h"
#include "prove.hpp"
#include "subcircuit.hpp"
#include "thread_pool.hpp"

#include <cadical.hpp>
#include <sstream>
#include <string>
#include <unordered_map>

using namespace talisman;

namespace talisman {

using entry_vec = std::vector<std::vector<std::pair<int32_t, int32_t>>>;
using thread_pool = std::unique_ptr<BS::thread_pool<BS::none>>;
extern thread_pool pool;
extern thread_local double local_guess_time;

// Collected SAT models for constant-column counterexamples (constant + inputs +
// nodes)
static std::vector<std::vector<int>> sat_repair_samples;

// For each constant column (excluding column 0) check via SAT whether the
// opposite Boolean value for the corresponding gate is reachable. On SAT,
// store the model in `sat_repair_samples` (do not modify the matrix in-place).
static void
check_constant_columns_sat(Subcircuit& sc,
                           nmod_mat_t M,
                           const std::string& label) {
  std::vector<int> const_cols = find_constant_columns(M);
  if(const_cols.empty()) {
    msg("%s: no constant columns found", label.c_str());
    return;
  }

  for(int col : const_cols) {
    if(col == 0)
      continue;// skip constant column

    Gate* g = nullptr;
    int n_inputs = (int)sc.inputs.size();
    if(col >= 1 && col <= n_inputs) {
      int idx = col - 1;
      int pos = (int)sc.inputs.size() - 1 - idx;
      if(pos >= 0 && pos < (int)sc.inputs.size())
        g = sc.inputs[pos];
    } else {
      int idx = col - 1 - n_inputs;
      int pos = (int)sc.nodes.size() - 1 - idx;
      if(pos >= 0 && pos < (int)sc.nodes.size())
        g = sc.nodes[pos];
    }

    if(!g) {
      warning("%s: constant column %d maps to no gate", label.c_str(), col);
      continue;
    }

    int base = nmod_mat_entry(M, 0, col) != 0;
    int var = (int)g->get_var() >> 1;
    int want_lit = base ? -var : var;

    CaDiCaL::Solver solver;
    for(auto& gate : sc.nodes) {
      const std::vector<int> cnf = gate->get_cnf_encoding();
      for(int lit : cnf)
        solver.add(lit);
    }
    solver.add(want_lit);
    solver.add(0);

    int res = solver.solve();
    if(res == 10) {
      msg("%s: column %d (gate %s) constant=%d, opposite reachable (SAT)",
          label.c_str(),
          col,
          g->get_var_name().c_str(),
          base);

      std::vector<int> model;
      model.reserve(get_size(sc) + 1);
      model.push_back(1);
      for(auto it = sc.inputs.rbegin(); it != sc.inputs.rend(); ++it) {
        int v = (*it)->get_var() >> 1;
        model.push_back(solver.val(v) > 0);
      }

      for(auto it = sc.nodes.rbegin(); it != sc.nodes.rend(); it++) {
        int var = (*it)->get_var() >> 1;
        int val = solver.val(var) > 0;
        model.push_back(val);
      }

      sat_repair_samples.push_back(std::move(model));
    } else if(res == 20) {
      msg("%s: column %d (gate %s) constant=%d, opposite NOT reachable (UNSAT)",
          label.c_str(),
          col,
          g->get_var_name().c_str(),
          base);
    } else {
      msg("%s: column %d (gate %s) constant=%d, solver returned %d",
          label.c_str(),
          col,
          g->get_var_name().c_str(),
          base,
          res);
    }
  }
}

/*------------------------------------------------------------------------*/

// Check a matrix (modular or packed bit-matrix) for constant columns
// and emit a formatted warning with the provided label.
static void
check_constant_columns(nmod_mat_t M, const std::string& label) {
  std::vector<int> const_cols = find_constant_columns(M);
  if(!const_cols.empty()) {
    std::ostringstream oss;
    for(size_t i = 0; i < const_cols.size(); ++i) {
      if(i)
        oss << ", ";
      oss << const_cols[i];
    }
    warning("%s: %s", label.c_str(), oss.str().c_str());
  } else {
    msg("%s: no constant columns found", label.c_str());
  }
}

/*------------------------------------------------------------------------*/

// Print an nmod matrix to the log with a label
static void
print_matrix(nmod_mat_t M, const std::string& label) {
  if(!M)
    return;
  msg("%s (size %d x %d)", label.c_str(), M->r, M->c);
  for(int i = 0; i < nmod_mat_nrows(M); ++i) {
    std::ostringstream oss;
    for(int j = 0; j < nmod_mat_ncols(M); ++j) {
      if(j)
        oss << ' ';
      oss << nmod_mat_entry(M, i, j);
    }
    msg("%s", oss.str().c_str());
  }
}
/*------------------------------------------------------------------------*/
// construct linear polynomials
bool
guesses_from_matrix_new(nmod_mat_t k, entry_vec& entries) {

  // check if we have polynomial for root gate
  if(k->r == 0 or nmod_mat_entry(k, k->r - 1, k->c - 1) == 0)
    return false;

  entries.resize(nmod_mat_nrows(k));

  for(int i = nmod_mat_nrows(k) - 1; i >= 0; --i) {
    auto row = nmod_mat_entry_ptr(k, i, 0);
    for(int j = nmod_mat_ncols(k) - 1; j >= 0; --j) {
      int32_t c = row[j];
      if(c != 0)
        entries[i].emplace_back(j, c);
    }
  }
  return true;
}
/*------------------------------------------------------------------------*/
// construct linear polynomials
bool
guesses_from_matrix(Subcircuit& sc,
                    nmod_mat_t k,
                    std::vector<LinPolynomial>& polies) {
  polies.clear();

  // check if we have polynomial for root gate
  if(k->r == 0 or nmod_mat_entry(k, k->r - 1, k->c - 1) == 0)
    return false;

  std::vector<int32_t> coeffs;
  std::vector<Gate*> mons;
  polies.reserve(nmod_mat_nrows(k));

  for(int i = nmod_mat_nrows(k) - 1; i >= 0; --i) {
    auto row = nmod_mat_entry_ptr(k, i, 0);
    for(int j = nmod_mat_ncols(k) - 1; j >= 0; --j) {
      int32_t c = row[j];
      if(c != 0) {
        Gate* g = get_gate(sc, j);
        c = 2 * c > k->mod.n ? c - k->mod.n : c;
        coeffs.push_back(c);
        mons.push_back(g);
      }
    }
    polies.emplace_back(coeffs, mons, true);
    coeffs.clear();
    mons.clear();
  }
  return true;
}

/*------------------------------------------------------------------------*/

void
complete_candidate(Subcircuit& sc, int32_t prime, LinPolynomial& poly) {

  nmod_mat_t k;
  nmod_mat_t v;
  nmod_mat_t kT;
  nmod_mat_t kTv;

  nmod_mat_init(v, 1, nmod_mat_nrows(k), prime);
  int i = 0;
  nmod_mat_set_entry(v, 0, i++, poly.get_constant());
  for(auto it = sc.inputs.rbegin(); it != sc.inputs.rend(); ++it)
    nmod_mat_set_entry(v, 0, i++, poly.get_coeff(*it));
  for(auto it = sc.nodes.rbegin(); it != sc.nodes.rend(); ++it) {
    nmod_mat_set_entry(v, 0, i++, poly.get_coeff(*it));
  }

  // normalize v to lc 1
  int c = 0;
  for(int j = v->c - 1; j >= 0; j--) {
    c = nmod_mat_entry(v, 0, j);
    if(c != 0)
      break;
  }
  c = nmod_inv(c, v->mod);
  for(int j = 0; j < v->c; j++) {
    int d = c * nmod_mat_entry(v, 0, j);
    nmod_mat_set_entry(v, 0, j, d);
  }

  nmod_mat_init(kT, k->c, k->r, prime);
  nmod_mat_init(kTv, k->c + 1, k->r, prime);
  nmod_mat_transpose(kT, k);
  nmod_mat_concat_vertical(kTv, kT, v);

  // mark pivots
  std::vector<int> pivots(kT->c, -1);
  for(int i = kT->r - 1; i >= 0; --i) {
    for(int j = kT->c - 1; j >= 0; j--) {
      if(nmod_mat_entry(kT, i, j) != 0) {
        pivots[j] = i;
        break;
      }
    }
  }

  // iterate over vector
  int I = kTv->r - 1;
  for(int j = kTv->c - 1; j >= 0; --j) {
    int c = nmod_mat_entry(kTv, I, j);
    if(c != 0 and pivots[j] >= 0)
      xmay(kTv, I, pivots[j], c);
  }

  // compute v - kTv[-1]
  for(int j = 0; j < kTv->c; j++) {
    int c = (int)nmod_mat_entry(v, 0, j) - (int)nmod_mat_entry(kTv, I, j);
    c = c % prime;
    c = c < 0 ? c + prime : c;
    nmod_mat_set_entry(kTv, I, j, c);
  }

  nmod_mat_clear(k);
  nmod_mat_init(k, kTv->c, kTv->r, prime);
  nmod_mat_transpose(k, kTv);

  nmod_mat_clear(kTv);
  nmod_mat_clear(kT);
  nmod_mat_clear(v);
}

/*------------------------------------------------------------------------*/

std::pair<bool, std::vector<LinPolynomial>>
guess_naive(Subcircuit& sc, nmod_mat_t mat_p) {

  nmod_mat_t ker;
  kernel(mat_p, ker);

  std::vector<LinPolynomial> guesses;
  bool candidate = guesses_from_matrix(sc, ker, guesses);

  nmod_mat_clear(ker);

  return std::make_pair(candidate, guesses);
}
/*------------------------------------------------------------------------*/
void
combine_multimod_polies(Subcircuit& sc,
                        std::vector<entry_vec>& entries,
                        std::vector<LinPolynomial>& polies) {
  polies.clear();
  polies.reserve(entries[0].size());

  for(int poly_idx = 0; poly_idx < entries[0].size(); poly_idx++) {
    LinPolynomial p;
    for(size_t i = 0; i < entries.size(); i++) {
      auto indices = entries[i][poly_idx];
      for(auto& [m, c] : indices) {
        Gate* g = get_gate(sc, m);
        p.add_term(g, c, i);
      }
    }
    polies.push_back(p);
  }
}
/*------------------------------------------------------------------------*/

bool
linear_polies_from_truth_table(Subcircuit& sc,
                               std::vector<LinPolynomial>& polies) {

  

  nmod_mat_t mat;
  sample_full_truth_table(sc, mat);
  if(verbose > 0)
    msg("Starting linear algebra for matrix of size %d x %d", mat->r, mat->c);

  std::vector<entry_vec> entries(PRIMES.size());
  for(size_t i = 0; i < PRIMES.size(); i++) {

    int32_t p = PRIMES[i];

    auto task = [&mat, &entries, p, i]() {
      // linear algebra
      nmod_mat_t ker_p;
      nmod_mat_t mat_p;

      nmod_mat_init_set(mat_p, mat);
      nmod_mat_set_mod(mat_p, p);

      kernel(mat_p, ker_p);

      // construct linear polynomials
      guesses_from_matrix_new(ker_p, entries[i]);

      nmod_mat_clear(mat_p);
      nmod_mat_clear(ker_p);
    };
    pool->detach_task(task);
  }
  pool->wait();

  nmod_mat_clear(mat);




  bool result; 
  if(std::any_of(
       entries.begin(), entries.end(), [](entry_vec& v) { return v.empty(); }))
    result = false;
  else {
      combine_multimod_polies(sc, entries, polies);
      result = true;
  }

  
  

  return result;
}
/*------------------------------------------------------------------------*/
thread_local std::vector<LinPolynomial> new_guesses;

std::pair<bool, std::vector<LinPolynomial>>
repair_guess(Subcircuit& sc,
             std::vector<std::vector<int>>& samples,
             nmod_mat_t mat) {

  if(verbose > 0)
    msg("Inserting %zu repair samples into mat", samples.size());

  append_repair_samples(mat, samples);

  // compute new kernel
  nmod_mat_t ker;
  kernel(mat, ker, samples.size());
  samples.clear();

  // construct linear polynomials
  new_guesses.clear();
  bool candidate = guesses_from_matrix(sc, ker, new_guesses);

  nmod_mat_clear(ker);

  return std::make_pair(candidate, new_guesses);
}

}// namespace talisman
