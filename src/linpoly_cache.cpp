#include "linpoly_cache.hpp"
#include "signal_statistics.hpp"

#include <unordered_map>
#include <vector>
#include <mutex>
#include <cadical.hpp>


namespace talisman {

namespace {
  std::unordered_map<const Gate*, std::vector<LinPolynomial>> g_map;
  const std::vector<LinPolynomial> kEmpty;
  std::mutex g_map_mutex;  // Mutex for thread-safe access
}

void linpoly_cache_init() {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  g_map.clear();
}

void linpoly_cache_clear() {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  g_map.clear();
}

void linpoly_cache_add(const Gate* g, const LinPolynomial& lp) {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  g_map[g].push_back(lp);
}

bool linpoly_cache_has(const Gate* g) {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  return g_map.find(g) != g_map.end();
}

const std::vector<LinPolynomial>& linpoly_cache_get(const Gate* g) {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  auto it = g_map.find(g);
  if (it == g_map.end()) return kEmpty;
  return it->second;
}


bool print_linpoly_cache() {
    if (verbose < 1) return false;  // Only print when verbose >= 1

    std::lock_guard<std::mutex> lock(g_map_mutex);
    if (g_map.empty()) {
      std::printf("[linpoly_cache] empty\n");
      return false;
    }

    std::printf("[linpoly_cache] entries = %zu\n", g_map.size());
    for (const auto& kv : g_map) {
      const Gate* g = kv.first;
      const auto& polys = kv.second;

      int gate_id = 0;
      const char* gate_name = "<null-gate>";
      if (g) {
        gate_id = static_cast<int>(g->get_var() / 2);
        gate_name = g->get_var_name().c_str();
      }

      std::printf("gate id=%d name=%s, polys=%zu\n",
                  gate_id, gate_name, polys.size());

      for (size_t i = 0; i < polys.size(); ++i) {
        std::printf("  [%zu] ", i);
        polys[i].print();
        // print() already prints newline when verbose >= 1
      }
    }
    return true;
  }

void linpoly_cache_for_each(const std::function<void(const Gate*, const std::vector<LinPolynomial>&)>& cb) {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  for (const auto &kv : g_map) {
    cb(kv.first, kv.second);
  }
}

// Helper function to check if two LinPolynomials are equal
static bool linpoly_equals(const LinPolynomial& a, const LinPolynomial& b) {
  auto terms_a = a.get_terms();
  auto terms_b = b.get_terms();

  if (terms_a.size() != terms_b.size()) return false;

  for (size_t i = 0; i < terms_a.size(); ++i) {
    if (terms_a[i].first != terms_b[i].first) return false;
    if (terms_a[i].second != terms_b[i].second) return false;
  }
  return true;
}

void linpoly_cache_deduplicate() {
  std::lock_guard<std::mutex> lock(g_map_mutex);
  size_t total_before = 0;
  size_t total_after = 0;

  for (auto& kv : g_map) {
    std::vector<LinPolynomial>& polys = kv.second;
    total_before += polys.size();

    if (polys.size() <= 1) {
      total_after += polys.size();
      continue;
    }

    // Keep only unique polynomials
    std::vector<LinPolynomial> unique_polys;
    for (const auto& poly : polys) {
      bool is_duplicate = false;
      for (const auto& existing : unique_polys) {
        if (linpoly_equals(poly, existing)) {
          is_duplicate = true;
          break;
        }
      }
      if (!is_duplicate) {
        unique_polys.push_back(poly);
      }
    }

    polys = std::move(unique_polys);
    total_after += polys.size();
  }

  std::printf("[linpoly_cache_deduplicate] %zu -> %zu polynomials (removed %zu duplicates)\n",
              total_before, total_after, total_before - total_after);
}

} // namespace talisman