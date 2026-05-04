// Polynomial.cpp

#include "polynomial.hpp"
#include "gate.hpp"
#include "signal_statistics.hpp"

#include <unordered_map>
#include <unordered_set>

namespace talisman {

/*---------------- Static Members ----------------*/

bool
Polynomial::is_zero() const {
  return multimod_terms.empty();
}

/*---------------- Monomial Implementation ----------------*/

// Helper: Compare two gates using the monomial ordering criterion
static bool
gate_less(Gate* g1, Gate* g2) {
  if(!g1 && !g2)
    return false;
  if(!g1)
    return false;
  if(!g2)
    return true;

  if(g1->get_cone_id() != g2->get_cone_id())
    return g1->get_cone_id() > g2->get_cone_id();// decreasing
  if(g1->get_level() != g2->get_level())
    return g1->get_level() > g2->get_level();// decreasing
  return g1->get_var_name() < g2->get_var_name();
}

Monomial::Monomial(std::vector<Gate*> v)
  : vars(std::move(v)) {
  sort_vars();
}

void
Monomial::sort_vars() {
  std::sort(vars.begin(), vars.end(), gate_less);
}

bool
Monomial::operator<(const Monomial& other) const {
  size_t n = std::min(vars.size(), other.vars.size());
  for(size_t i = 0; i < n; ++i) {
    Gate* g1 = vars[i];
    Gate* g2 = other.vars[i];

    if(gate_less(g1, g2))
      return true;
    if(gate_less(g2, g1))
      return false;
  }
  // If one monomial is prefix of the other, the **longer monomial comes
  // first**.
  return vars.size() > other.vars.size();
}

void
Monomial::print() const {
  if(vars.empty()) {
    std::cout << "1" << std::endl;// constant term
    return;
  }

  for(size_t i = 0; i < vars.size(); ++i) {
    if(i > 0)
      std::cout << "*";
    if(vars[i]) {
      std::cout << vars[i]->get_var_name();
    } else {
      std::cout << "null";
    }
  }
  std::cout << std::endl;
}

Monomial
Monomial::get_remainder() const {
  // Return a monomial with all variables except the leading one
  if(vars.empty()) {
    return Monomial();// empty monomial
  }
  // Skip the first variable (the leading one)
  std::vector<Gate*> remainder_vars(vars.begin() + 1, vars.end());
  return Monomial(remainder_vars);
}

/*---------------- Polynomial Implementation ----------------*/

Polynomial::Polynomial(Gate* gate) {
  if(gate) {
    Monomial m({ gate });
    CoeffVec coeff;
    coeff.set_constant(1);
    multimod_terms[m] = coeff;
  }
}

Polynomial::Polynomial(int coeff) {
  if(coeff != 0) {
    Monomial m(std::vector<Gate*>{});// empty monomial = constant term
    CoeffVec cv;
    cv.set_constant(coeff);
    multimod_terms[m] = cv;
  }
}

/*-------------------------------------------------------------------------------*/
void
Polynomial::mul_inplace(const Polynomial& other) {
  if(multimod_terms.empty() || other.multimod_terms.empty()) {
    multimod_terms.clear();
    return;
  }

  MultimodTermMap result;
  std::vector<Gate*> merged;
  for(const auto& [a_m, a_c] : multimod_terms) {
    for(const auto& [b_m, b_c] : other.multimod_terms) {
      // Merge variable sets with single allocation
      merged.clear();
      merged.reserve(a_m.vars.size() + b_m.vars.size());
      merged.insert(merged.end(), a_m.vars.begin(), a_m.vars.end());
      merged.insert(merged.end(), b_m.vars.begin(), b_m.vars.end());

      // Remove duplicates (only one pass)
      std::sort(merged.begin(), merged.end(), [](Gate* g1, Gate* g2) {
        return g1 < g2;
      });
      merged.erase(std::unique(merged.begin(),
                               merged.end(),
                               [](Gate* g1, Gate* g2) { return g1 == g2; }),
                   merged.end());
      // Do NOT sort again - merged is already in pointer order

      Monomial m(std::move(merged));// Move constructor avoids copy
      CoeffVec prod = a_c * b_c;

      if(!prod.is_zero()) {
        auto [it, inserted] = result.try_emplace(m, prod);
        if(!inserted) {
          // Use in-place addition instead of element-wise loop
          it->second.add_inplace(prod);
          if(it->second.is_zero()) {
            result.erase(it);
          }
        }
      }
    }
  }

  multimod_terms.swap(result);
}

/*-------------------------------------------------------------------------------*/
std::string
Polynomial::to_string() const {
  if(multimod_terms.empty())
    return "0";

  std::ostringstream oss;
  bool first = true;

  for(const auto& [m, coeff] : multimod_terms) {
    if(coeff.is_zero())
      continue;

    if(!first) {
      oss << " + ";
    }

    // Print coefficient vector: {c1, c2, ...}
    oss << "{";
    for(size_t i = 0; i < coeff.entries.size(); ++i) {
      oss << coeff.entries[i];
      if(i + 1 < coeff.entries.size())
        oss << ", ";
    }
    oss << "}";

    // Print monomial variables
    if(!m.vars.empty()) {
      oss << "*";
      for(size_t i = 0; i < m.vars.size(); ++i) {
        if(i > 0)
          oss << "*";
        if(m.vars[i]) {
          oss << m.vars[i]->get_var_name();
        } else {
          oss << "null";
        }
      }
    }

    first = false;
  }

  if(first)
    return "0";// nothing non-zero printed
  return oss.str();
}
/*---------------------------------------------------------------------*/
void
Polynomial::print() const {
  std::cout << to_string() << std::endl;
}

/*---------------------------------------------------------------------*/
void
Polynomial::print_multimod() const {
  if(multimod_terms.empty()) {
    std::cout << "0" << std::endl;
    return;
  }

  bool first = true;
  for(const auto& [m, coeff] : multimod_terms) {
    if(coeff.is_zero())
      continue;

    if(!first) {
      std::cout << " + ";
    }

    coeff.print();

    if(!m.vars.empty()) {
      std::cout << "*";
      for(size_t i = 0; i < m.vars.size(); ++i) {
        if(i > 0)
          std::cout << "*";
        if(m.vars[i]) {
          std::cout << m.vars[i]->get_var_name();
        } else {
          std::cout << "null";
        }
      }
    }

    first = false;
  }

  std::cout << std::endl;
}

/*-------------------------------------------------------------------------------*/
bool
Polynomial::check_leading_coeff() {
  extern std::vector<CoeffVec> two_pow_NN_imgs;

  auto lc = multimod_terms.begin()->second;

  if(std::any_of(two_pow_NN_imgs.begin(),
                 two_pow_NN_imgs.end(),
                 [&lc](const CoeffVec& v) { return lc == v; })) {

    multimod_terms.erase(multimod_terms.begin());
    return true;
  }
  return false;
}

/*-------------------------------------------------------------------------------*/

// Remove or reject coefficients equal to 2^NN modulo all primes
void
Polynomial::check_coeffs(bool remove) {
  // vector is defined in linpolynomial.cpp
  extern std::vector<CoeffVec> two_pow_NN_imgs;

  if(!remove) {
    for(const auto& [mono, coeff] : multimod_terms) {
      if(std::any_of(two_pow_NN_imgs.begin(),
                     two_pow_NN_imgs.end(),
                     [&coeff](const CoeffVec& v) { return coeff == v; })) {
        throw std::runtime_error("coeff equals 2^NN modulo all primes");
      }
    }
    return;
  }

  std::vector<Monomial> to_erase;
  for(const auto& [mono, coeff] : multimod_terms) {
    if(std::any_of(two_pow_NN_imgs.begin(),
                   two_pow_NN_imgs.end(),
                   [&coeff](const CoeffVec& v) { return coeff == v; })) {
      to_erase.push_back(mono);
    }
  }

  for(const auto& mono : to_erase) {
    std::string name;
    if(mono.vars.empty()) {
      name = "constant term";
    } else {
      for(size_t i = 0; i < mono.vars.size(); ++i) {
        if(i)
          name += "*";
        name += mono.vars[i]->get_var_name();
      }
    }
    std::cerr << "Removing monomial " << name
              << " because coeff is multiple of 2^NN" << std::endl;
    multimod_terms.erase(mono);
  }
}

/*-------------------------------------------------------------------------------*/

// Scalar operators needed for gate polynomial generation
Polynomial
operator*(Polynomial a, const Polynomial& b) {
  a.mul_inplace(b);
  return a;
}

/*-------------------------------------------------------------------------------*/

Polynomial
operator-(Polynomial a, const Polynomial& b) {
  for(const auto& [m, c] : b.multimod_terms) {
    auto it = a.multimod_terms.find(m);
    if(it != a.multimod_terms.end()) {
      // Subtract coefficients element-wise
      for(size_t i = 0; i < it->second.entries.size(); ++i) {
        it->second.entries[i]
          = nmod_sub(it->second.entries[i], c.entries[i], moduli[i]);
      }
      if(it->second.is_zero()) {
        a.multimod_terms.erase(it);
      }
    } else {
      // Negate and add
      CoeffVec neg_c;
      for(size_t i = 0; i < c.entries.size(); ++i) {
        neg_c.entries[i] = nmod_neg(c.entries[i], moduli[i]);
      }
      a.multimod_terms[m] = neg_c;
    }
  }
  return a;
}

/*---------------- Gate Polynomial Generation ----------------*/

// Helper: polynomial from an AIG literal (always multimodular)
Polynomial
Polynomial::literal_to_poly(unsigned lit) {
  Polynomial result;

  CoeffVec coeff_one;
  coeff_one.set_constant(1);
  CoeffVec coeff_zero;
  coeff_zero.set_constant(0);

  // Constants
  if(lit == 0u) {
    // false = 0
    Monomial m(std::vector<Gate*>{});
    result.multimod_terms[m] = coeff_zero;
    return result;
  }
  if(lit == 1u) {
    // true = 1
    Monomial m(std::vector<Gate*>{});
    result.multimod_terms[m] = coeff_one;
    return result;
  }

  // Extract gate by node id
  Gate* child = getGateByNodeId(static_cast<int>(lit / 2));
  Monomial m({ child });
  result.multimod_terms[m] = coeff_one;

  const bool complemented = (lit & 1u) != 0u;
  if(complemented) {
    // (1 - x): constant term 1, variable term -1
    Monomial const_m(std::vector<Gate*>{});
    result.multimod_terms[const_m] = coeff_one;

    CoeffVec coeff_neg_one;
    coeff_neg_one.set_constant(-1);
    result.multimod_terms[m] = coeff_neg_one;
  }

  return result;
}

Polynomial
Polynomial::get_gate_poly(Gate* g) {
  if(!g)
    throw std::runtime_error("get_gate_poly: null gate");

  // Create multimodular gate polynomial
  Polynomial gout;
  CoeffVec coeff_one;
  coeff_one.set_constant(1);
  Monomial m({ g });
  gout.multimod_terms[m] = coeff_one;

  // Child polynomials respecting signedness/complementation
  Polynomial c0 = literal_to_poly(g->get_ch0());
  Polynomial c1 = literal_to_poly(g->get_ch1());

  // Nonlinear gate polynomial: g - c0*c1
  Polynomial prod = c0 * c1;
  Polynomial result = gout - prod;
  return result;
}

/*---------------- Polynomial Reduction ----------------*/

void
Polynomial::reduce_by(const Polynomial& p2) {
  if(!p2.multimod_terms.empty()) {
    return reduce_by_multimod(p2);
  } else {
    die(1, "Polynomial::reduce_by: p2 is empty");
  }
}

/*---------------- Multimodular Operations ----------------*/

void
Polynomial::add_term(const Monomial& mono, int32_t coeff, int idx) {
  if(coeff == 0)
    return;

  multimod_terms[mono].set_at_idx(coeff, idx);
}

void
Polynomial::add_term(const Monomial& mono, CoeffVec& coeff) {
  if(coeff.is_zero())
    return;

  auto it = multimod_terms.find(mono);

  if(it == multimod_terms.end()) {
    multimod_terms.emplace(mono, std::move(coeff));
  }
}

void
Polynomial::addmul_term(const Monomial& mono,
                        const CoeffVec& multiplier,
                        const CoeffVec& coeff) {
  if(coeff.is_zero())
    return;

  auto it = multimod_terms.find(mono);
  if(it != multimod_terms.end()) {
    CoeffVec& my_coeff = it->second;
    my_coeff.addmul(multiplier, coeff);

    if(my_coeff.is_zero()) {
      multimod_terms.erase(it);
    }
  } else {
    multimod_terms.emplace(mono, multiplier * coeff);
  }
}

std::vector<Gate*> combined_vars;
void
Polynomial::addmul_term(const Monomial& mono,
                        const Monomial& remainder,
                        const CoeffVec& multiplier,
                        const CoeffVec& coeff) {
  if(coeff.is_zero())
    return;

  // Merge variable sets efficiently
  combined_vars.clear();
  combined_vars.reserve(mono.vars.size() + remainder.vars.size());
  combined_vars.insert(combined_vars.end(), mono.vars.begin(), mono.vars.end());
  combined_vars.insert(
    combined_vars.end(), remainder.vars.begin(), remainder.vars.end());

  // Remove duplicates (single pass, maintain pointer order)
  std::sort(combined_vars.begin(), combined_vars.end(), [](Gate* g1, Gate* g2) {
    return g1 < g2;
  });
  combined_vars.erase(std::unique(combined_vars.begin(),
                                  combined_vars.end(),
                                  [](Gate* g1, Gate* g2) { return g1 == g2; }),
                      combined_vars.end());
  // Pointer-sorted is already gate_less-compatible, no need for sort_vars()

  Monomial combined_mono(std::move(combined_vars));
  CoeffVec prod = multiplier * coeff;

  auto [it, inserted] = multimod_terms.try_emplace(combined_mono, prod);
  if(!inserted) {
    it->second.add_inplace(prod);
    if(it->second.is_zero()) {
      multimod_terms.erase(it);
    }
  }
}

std::vector<std::pair<Monomial, CoeffVec>> reducible_terms;
std::vector<Gate*> remainder_vars;
void
Polynomial::reduce_by_multimod(const Polynomial& p2) {
  const bool do_verbose = verbose > 1;

  if(p2.multimod_terms.empty())
    throw std::runtime_error("Cannot reduce by empty multimod polynomial");

  const auto other_lead = *p2.multimod_terms.begin();
  if(other_lead.first.vars.empty() || other_lead.first.vars[0] == nullptr)
    throw std::runtime_error(
      "Leading monomial of reducer must have a variable");

  nonlin_steps++;

  Gate* target_var = other_lead.first.vars[0];

  CoeffVec other_coeff = other_lead.second;

  // Collect all leading terms that still contain the target variable; map is
  // ordered
  reducible_terms.clear();
  reducible_terms.reserve(multimod_terms.size());
  for(const auto& kv : multimod_terms) {
    const Monomial& mono = kv.first;
    bool contains_target = false;
    for(auto* g : mono.vars) {
      if(g == target_var) {
        contains_target = true;
        break;
      }
    }

    if(!contains_target)
      break;// abort search once a term without target var is seen
    reducible_terms.push_back(kv);
  }

  if(reducible_terms.empty()) {
    if(do_verbose)
      std::cout << "reduce_by_multimod: no term contains target var\n";
    return;
  }

  if(do_verbose) {
    std::cout << "\n=== reduce_by_multimod START ===" << std::endl;
    std::cout << "Reducing (multimod): " << this->to_string() << std::endl;
    std::cout << "By p2 (multimod):    " << p2.to_string() << std::endl;
    std::cout << "Other leading coeff: ";
    other_coeff.print();
    std::cout << std::endl;
  }

  for(const auto& [mono, my_coeff] : reducible_terms) {
    if(my_coeff.is_zero())
      continue;

    // Build remainder by removing the target variable from the monomial
    remainder_vars.clear();
    remainder_vars.reserve(mono.vars.size());
    bool removed = false;
    for(auto* g : mono.vars) {
      if(!removed && g == target_var) {
        removed = true;
        continue;
      }
      remainder_vars.push_back(g);
    }
    Monomial remainder(std::move(remainder_vars));

    CoeffVec multiplier = my_coeff.negdiv(other_coeff);

    if(do_verbose) {
      std::cout << "Reducing term:       ";
      my_coeff.print();
      std::cout << " * ";
      mono.print();
      std::cout << "Multiplier:          ";
      multiplier.print();
      std::cout << "Remainder:           ";
      remainder.print();
    }

    for(const auto& [p2_mono, p2_coeff] : p2.multimod_terms) {
      addmul_term(p2_mono, remainder, multiplier, p2_coeff);
    }
  }

  // Drop any terms whose coefficients equal 2^NN across all primes
  // check_coeffs(true);

  if(do_verbose)
    std::cout << "=== reduce_by_multimod END ===" << std::endl << std::endl;
}

/*---------------- LinPolynomial Conversion ----------------*/

Polynomial
Polynomial::from_linpoly(const LinPolynomial& linpoly) {
  Polynomial result;

  if(linpoly.multimod) {
    // Pre-create constant monomial to avoid repeated allocations
    Monomial constant_monomial(std::vector<Gate*>{});

    for(const auto& [gate, coeff_vec] : linpoly.multimod_terms) {
      if(gate == nullptr) {
        result.multimod_terms[constant_monomial] = coeff_vec;
      } else {
        Monomial linear_monomial({ gate });
        result.multimod_terms[linear_monomial] = coeff_vec;
      }
    }
  } else {
    // Convert single-modulus terms to multimod
    auto terms_list = linpoly.get_terms();

    for(const auto& [gate, coeff] : terms_list) {
      if(coeff == 0)
        continue;

      // Convert single coefficient to multimodular coefficient vector
      CoeffVec coeff_vec;
      coeff_vec.set_constant(coeff);

      if(gate == nullptr) {
        Monomial constant_monomial(std::vector<Gate*>{});
        result.multimod_terms[constant_monomial] = coeff_vec;
      } else {
        Monomial linear_monomial({ gate });
        result.multimod_terms[linear_monomial] = coeff_vec;
      }
    }
  }

  return result;
}

std::unordered_set<Gate*> variables;
std::vector<std::pair<Monomial, CoeffVec>> terms_to_add;
std::vector<Monomial> terms_to_remove;
std::vector<Gate*> new_vars;

void
Polynomial::reduce_single_children() {
  if(verbose >= 2) {
    msg("=== Starting reduce_single_children ===");
  }

  bool reduced = true;
  int counter = 0;

  while(reduced) {
    reduced = false;

    if(multimod_terms.empty())
      break;

    // Get the leading monomial to exclude it from the check
    const Monomial& leading_mono = multimod_terms.begin()->first;

    // Collect all variables (gates) in all monomials, excluding leading term
    // Use unordered_set instead of std::set for O(1) lookups
    variables.clear();
    for(const auto& [mono, coeff] : multimod_terms) {
      // Cache comparison result to avoid double comparison
      bool is_leading = !(mono < leading_mono) && !(leading_mono < mono);
      if(!is_leading) {
        for(Gate* g : mono.vars) {
          if(g)
            variables.insert(g);
        }
      }
    }

    if(verbose >= 2) {
      msg("Checking %zu variables for single-parent reduction (excluding "
          "leading term)",
          variables.size());
    }

    // Check each variable (gate) to see if it has only one parent
    for(Gate* child_gate : variables) {
      if(!child_gate)
        continue;

      size_t child_parent_size = child_gate->parent_size();

      // If child has exactly one parent, reduce by child's gate polynomial
      if(child_parent_size == 1) {
        if(verbose >= 2) {
          msg("Gate %s has exactly 1 parent, reducing by its gate polynomial",
              child_gate->get_var_name().c_str());
        }

        // Get the gate polynomial for this child
        Polynomial child_poly = Polynomial::get_gate_poly(child_gate);

        if(verbose >= 2) {
          msg("Reducing by gate polynomial of %s:",
              child_gate->get_var_name().c_str());
          child_poly.print_multimod();
        }

        // Custom reduction: reduce all terms that contain child_gate
        if(child_poly.multimod_terms.empty()) {
          if(verbose >= 2) {
            msg("Child polynomial is empty, skipping reduction");
          }
          continue;
        }

        const auto& child_leading = *child_poly.multimod_terms.begin();
        Gate* child_leading_var = child_leading.first.vars.empty()
                                    ? nullptr
                                    : child_leading.first.vars[0];
        const CoeffVec& child_leading_coeff = child_leading.second;

        if(!child_leading_var) {
          if(verbose >= 2) {
            msg("Child polynomial has no leading variable, skipping");
          }
          continue;
        }

        // Batch collect terms to modify
        terms_to_add.clear();
        terms_to_remove.clear();
        terms_to_add.reserve(multimod_terms.size());

        for(const auto& [mono, coeff] : multimod_terms) {
          // Find position of child variable
          auto it
            = std::find(mono.vars.begin(), mono.vars.end(), child_leading_var);
          if(it == mono.vars.end())
            continue;// Not in this monomial

          size_t child_var_pos = std::distance(mono.vars.begin(), it);

          if(verbose >= 2) {
            msg("Found term containing %s, reducing it",
                child_leading_var->get_var_name().c_str());
          }

          // Build remainder monomial efficiently
          remainder_vars.clear();
          remainder_vars.reserve(mono.vars.size() - 1);
          for(size_t i = 0; i < mono.vars.size(); ++i) {
            if(i != child_var_pos) {
              remainder_vars.push_back(mono.vars[i]);
            }
          }
          Monomial remainder(std::move(remainder_vars));

          // Calculate multiplier: -coeff / child_leading_coeff
          CoeffVec multiplier = coeff.negdiv(child_leading_coeff);

          // Add multiplied terms from child_poly to the result
          for(const auto& [child_mono, child_coeff] :
              child_poly.multimod_terms) {
            // Cache leading check
            bool is_leading = !(child_mono < child_leading.first)
                              && !(child_leading.first < child_mono);
            if(is_leading)
              continue;

            // Build new monomial efficiently
            new_vars.clear();
            new_vars.reserve(remainder.vars.size() + child_mono.vars.size());
            new_vars.insert(
              new_vars.end(), remainder.vars.begin(), remainder.vars.end());
            new_vars.insert(
              new_vars.end(), child_mono.vars.begin(), child_mono.vars.end());

            // Remove duplicates (single pass, pointer order maintained)
            std::sort(new_vars.begin(), new_vars.end(), [](Gate* g1, Gate* g2) {
              return g1 < g2;
            });
            new_vars.erase(
              std::unique(new_vars.begin(),
                          new_vars.end(),
                          [](Gate* g1, Gate* g2) { return g1 == g2; }),
              new_vars.end());

            Monomial new_mono(std::move(new_vars));

            // Multiply: multiplier * child_coeff
            CoeffVec new_coeff = multiplier * child_coeff;
            terms_to_add.push_back({ new_mono, new_coeff });
          }

          // Mark original term for removal
          terms_to_remove.push_back(mono);
          reduced = true;
        }

        // Remove the reduced terms
        for(const auto& mono : terms_to_remove) {
          multimod_terms.erase(mono);
        }

        // Add the new terms using batch insertion
        for(const auto& [mono, coeff] : terms_to_add) {
          if(!coeff.is_zero()) {
            auto [it, inserted] = multimod_terms.try_emplace(mono, coeff);
            if(!inserted) {
              it->second.add_inplace(coeff);
              if(it->second.is_zero()) {
                multimod_terms.erase(it);
              }
            }
          }
        }

        if(verbose >= 2) {
          msg("Size after reduction: %zu", this->get_size());
          msg("Polynomial after reducing by %s:",
              child_gate->get_var_name().c_str());
          this->print_multimod();
        }

        if(reduced) {
          break;// Restart the loop since polynomial structure may have changed
        }
      }
    }
  }

  if(verbose >= 2) {
    msg("=== Completed reduce_single_children ===");
    msg("Final size: %zu, counter: %d", this->get_size(), counter);
  }
}

void
translate_all_stored_gate_polys() {
  for(int i = 0; i < (int)num_gates; i++) {
    Gate* g = gates[i];
    if(g->has_gatepoly()) {
      LinPolynomial* linpoly = g->get_gatepoly_mut();
      Polynomial poly = Polynomial::from_linpoly(*linpoly);

      // check whether leading term of poly matches leading_gate
      auto it = poly.multimod_terms.begin();
      if(it == poly.multimod_terms.end()) {
        throw std::runtime_error("Stored polynomial is empty");
      }
      Monomial leading_mono = it->first;
      // store the polynomial to the proper gate. This can also be now a
      // different gate if translation changed the leading variable
      if(leading_mono.vars.empty() || leading_mono.vars[0] == nullptr) {
        throw std::runtime_error(
          "Leading monomial of stored polynomial must have a variable");
      }

      Gate* poly_leading_gate = leading_mono.vars[0];

      // If the leading variable changed due to translation, move polynomial to
      // the new leading gate
      if(poly_leading_gate != g) {
        if(verbose > 0)
          msg("Translation changed leading variable from %s to %s",
              g->get_var_name().c_str(),
              poly_leading_gate->get_var_name().c_str());

        // Move the polynomial to the new leading gate
        poly_leading_gate->set_gatepoly(*linpoly);

        // Clear from the old gate
        g->clear_gatepoly();
      }
    }
  }
}

}// namespace talisman
