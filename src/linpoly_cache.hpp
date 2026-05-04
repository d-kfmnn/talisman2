#ifndef TALISMAN_SRC_LINPOLY_CACHE_HPP_
#define TALISMAN_SRC_LINPOLY_CACHE_HPP_

#include <vector>
#include <unordered_map>
#include <cstdint>

#include "linpolynomial.hpp" 
#include "gate.hpp"
#include <functional>

namespace talisman {

void linpoly_cache_init();

void linpoly_cache_clear();

void linpoly_cache_add(const Gate* g, const LinPolynomial& lp);

bool linpoly_cache_has(const Gate* g);

const std::vector<LinPolynomial>& linpoly_cache_get(const Gate* g);

bool print_linpoly_cache();

// Iterate all entries in the linpoly cache. The callback receives (Gate*, polys).
void linpoly_cache_for_each(const std::function<void(const Gate*, const std::vector<LinPolynomial>&)>& cb);

// Remove duplicate polynomials from the cache (for each gate, keep only unique polys)
void linpoly_cache_deduplicate();

} // namespace talisman

#endif // TALISMAN_SRC_LINPOLY_CACHE_HPP_