/*------------------------------------------------------------------------*/
/*! \file prove.hpp
    \brief contains the functions to prove the guess

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_PROVE_H_
#define TALISMAN_SRC_PROVE_H_
/*------------------------------------------------------------------------*/
#include "linpolynomial.hpp"
#include "subcircuit.hpp"
#include <cadical.hpp>
#include <chrono>
#include <list>

namespace talisman {
extern std::chrono::duration<double> total_sat_time;
extern std::chrono::duration<double> total_enc_poly_time;
extern std::chrono::duration<double> total_enc_circuit_time;
extern std::chrono::duration<double> total_prove_time;

int32_t
encode_subcircuit(Subcircuit& sc, CaDiCaL::Solver* solver);

bool
prove_guess(Subcircuit& sc,
            LinPolynomial& guess,
            CaDiCaL::Solver* solver,
            std::vector<std::vector<int>>& repair_samples,
            int prime,
            int32_t& fresh_var,
            std::vector<int> factors = { 0 });

bool
check_linpoly_cache(int sc_depth,
                    int prime,
                    std::vector<int> factors = { 0, 1 });

// When trust_abc is enabled, trust the linear polynomial cache produced by
// ABC and attach cached polynomials directly to gates without SAT checks.
bool
trust_linpoly_cache();

#ifdef USE_GUROBI

void
start_MIP();

#endif

}// namespace talisman

#endif// TALISMAN_SRC_PROVE_H_
