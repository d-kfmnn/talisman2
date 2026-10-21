/*------------------------------------------------------------------------*/
/*! \file candidates.cpp
    \brief contains functions necessary to detect half-/full-adders and the
           reverse "free simplification area" boundary directly on the
           parsed AIG

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "candidates.hpp"

#include "aig.hpp"
#include "signal_statistics.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/*------------------------------------------------------------------------*/
// This file replaces the "supertiny fraction" of ABC that TalisMan used to
// depend on for candidate detection:
//   - structural half-/full-adder recognition (ABC's acecRe.c/acecXor.c),
//   - the reverse local-K-cover boundary used to bound guessing/repair.
//
// Both are re-derived directly from the AIG already parsed into aig.cpp's
// 'model' by parse_aig(), so no second parse and no external EDA library are
// needed. The algorithms mirror ABC's (bounded k<=3 cut enumeration paired
// by shared leaves for adders; PO-driven BFS/DFS for the reverse cover) but
// classify cuts via direct truth-table simulation instead of ABC's packed
// truth-table tables, which is both easier to verify and avoids ABC's
// generic Vec_Int_t/hash-table machinery.
/*------------------------------------------------------------------------*/

namespace talisman {
namespace {

/*------------------------------------------------------------------------*/
// Lightweight AIG view, built once directly from aig.cpp's parsed model.
/*------------------------------------------------------------------------*/
struct AigView {
  unsigned M = 0;// number of variables; ids are 0..M-1, 0 = constant
  std::vector<int> fanin0, fanin1;
  std::vector<uint8_t> compl0, compl1;
  std::vector<uint8_t> is_and;
  std::vector<int> refcount;// fanout count, used by the trivial-HA filter
  std::vector<int> level;   // logic level, used by the K-cover boundary
  bool built = false;
};

AigView g_view;

void
ensure_aig_view() {
  if(g_view.built)
    return;
  aiger* model = get_model();
  unsigned M = model->maxvar + 1;
  g_view.M = M;
  g_view.fanin0.assign(M, 0);
  g_view.fanin1.assign(M, 0);
  g_view.compl0.assign(M, 0);
  g_view.compl1.assign(M, 0);
  g_view.is_and.assign(M, 0);
  g_view.refcount.assign(M, 0);
  g_view.level.assign(M, 0);

  // The AIGER binary format guarantees 'ands' is stored in topological
  // order (every rhs var is smaller than lhs var), so a single forward pass
  // is enough to fill fanins and levels.
  for(unsigned i = 0; i < model->num_ands; i++) {
    const aiger_and& a = model->ands[i];
    unsigned v = a.lhs >> 1;
    if(v >= M)
      continue;
    unsigned f0 = a.rhs0 >> 1, f1 = a.rhs1 >> 1;
    g_view.is_and[v] = 1;
    g_view.fanin0[v] = static_cast<int>(f0);
    g_view.compl0[v] = static_cast<uint8_t>(a.rhs0 & 1u);
    g_view.fanin1[v] = static_cast<int>(f1);
    g_view.compl1[v] = static_cast<uint8_t>(a.rhs1 & 1u);
    g_view.refcount[f0]++;
    g_view.refcount[f1]++;
    int lv0 = g_view.level[f0], lv1 = g_view.level[f1];
    g_view.level[v] = 1 + (lv0 > lv1 ? lv0 : lv1);
  }
  for(unsigned i = 0; i < model->num_outputs; i++) {
    unsigned v = model->outputs[i].lit >> 1;
    if(v < M)
      g_view.refcount[v]++;
  }
  g_view.built = true;
}

/*------------------------------------------------------------------------*/
// k<=3-leaf cuts with exact truth tables (no ABC-style packed byte tricks).
/*------------------------------------------------------------------------*/
struct Cut {
  uint8_t n = 0;// number of leaves: 1, 2 or 3
  int leaves[3] = { 0, 0, 0 };// sorted ascending
  unsigned truth = 0;// truth table over the low 2^n bits
};

// Truth table (over 2^n rows) of "variable i" among n ordered variables.
unsigned
base_pattern(int i, int n) {
  unsigned t = 0;
  for(int row = 0; row < (1 << n); row++)
    if((row >> i) & 1)
      t |= (1u << row);
  return t;
}

// Re-expresses cut c's truth table in terms of a (super-)set of leaves.
unsigned
expand_truth(const Cut& c, const int* new_leaves, int new_n) {
  unsigned result = 0;
  for(int row = 0; row < (1 << new_n); row++) {
    int old_row = 0;
    for(int i = 0; i < c.n; i++)
      for(int j = 0; j < new_n; j++)
        if(new_leaves[j] == c.leaves[i]) {
          if((row >> j) & 1)
            old_row |= (1 << i);
          break;
        }
    if((c.truth >> old_row) & 1u)
      result |= (1u << row);
  }
  return result;
}

// Merges two (<=3-leaf) cuts into a new cut for the AND node that combines
// them; fails (returns false) if the union of leaves would exceed 3.
bool
merge_cuts(const Cut& c0, bool compl0, const Cut& c1, bool compl1, Cut& out) {
  int merged[3];
  int m = 0;
  for(int i = 0; i < c0.n; i++)
    merged[m++] = c0.leaves[i];
  for(int i = 0; i < c1.n; i++) {
    bool found = false;
    for(int j = 0; j < m; j++)
      if(merged[j] == c1.leaves[i]) {
        found = true;
        break;
      }
    if(!found) {
      if(m == 3)
        return false;
      merged[m++] = c1.leaves[i];
    }
  }
  std::sort(merged, merged + m);
  out.n = static_cast<uint8_t>(m);
  for(int i = 0; i < m; i++)
    out.leaves[i] = merged[i];

  unsigned mask = (1u << (1 << m)) - 1u;
  unsigned t0 = expand_truth(c0, out.leaves, m);
  unsigned t1 = expand_truth(c1, out.leaves, m);
  if(compl0)
    t0 = (~t0) & mask;
  if(compl1)
    t1 = (~t1) & mask;
  out.truth = t0 & t1 & mask;
  return true;
}

bool
cut_equal(const Cut& a, const Cut& b) {
  if(a.n != b.n || a.truth != b.truth)
    return false;
  for(int i = 0; i < a.n; i++)
    if(a.leaves[i] != b.leaves[i])
      return false;
  return true;
}

/*------------------------------------------------------------------------*/
// Phase tables: which (input polarity, output polarity) explains a given
// raw AND2 (n=2) / MAJ3 (n=3) truth table. Built once, generated directly
// from the definitions rather than hand-transcribed constants.
/*------------------------------------------------------------------------*/
struct PhaseInfo {
  int input_phase[3] = { 0, 0, 0 };
  int output_phase = 0;
};

std::unordered_map<unsigned, PhaseInfo>
build_maj_table(int n) {
  std::unordered_map<unsigned, PhaseInfo> table;
  unsigned mask = (1u << (1 << n)) - 1u;
  unsigned b[3];
  for(int i = 0; i < n; i++)
    b[i] = base_pattern(i, n);

  for(int phase_mask = 0; phase_mask < (1 << n); phase_mask++) {
    unsigned p[3];
    for(int i = 0; i < n; i++)
      p[i] = ((phase_mask >> i) & 1) ? ((~b[i]) & mask) : b[i];
    unsigned fn = (n == 2) ? (p[0] & p[1])
                           : ((p[0] & p[1]) | (p[1] & p[2]) | (p[0] & p[2]));
    for(int out_c = 0; out_c < 2; out_c++) {
      unsigned t = out_c ? ((~fn) & mask) : fn;
      // Majority is self-dual (flipping every input complements the
      // output), so a handful of (phase, out_c) pairs alias to the same
      // truth value; keeping the first one found is an arbitrary but
      // consistent and equally valid choice.
      PhaseInfo info;
      for(int i = 0; i < n; i++)
        info.input_phase[i] = (phase_mask >> i) & 1;
      info.output_phase = out_c;
      table.emplace(t, info);
    }
  }
  return table;
}

const std::unordered_map<unsigned, PhaseInfo>&
maj_table(int n) {
  static const auto t2 = build_maj_table(2);
  static const auto t3 = build_maj_table(3);
  return n == 2 ? t2 : t3;
}

unsigned
xor_pattern(int n, const int* input_phase) {
  unsigned mask = (1u << (1 << n)) - 1u;
  unsigned t = 0;
  for(int i = 0; i < n; i++) {
    unsigned b = base_pattern(i, n);
    if(input_phase[i])
      b = (~b) & mask;
    t ^= b;
  }
  return t & mask;
}

/*------------------------------------------------------------------------*/
// Candidate half-/full-adder pairs: a k-leaf XOR cut and a k-leaf AND/MAJ
// cut that share exactly the same leaves.
/*------------------------------------------------------------------------*/
struct AdderCandidate {
  int n_in = 0;// 2 (half-adder) or 3 (full-adder)
  int in[3] = { 0, 0, 0 };
  int xor_root = 0, maj_root = 0;
  int input_phase[3] = { 0, 0, 0 };
  int output_phase_xor = 0, output_phase_maj = 0;
  unsigned xor_truth = 0, maj_truth = 0;
};

struct LeafKey {
  int a = 0, b = 0, c = 0;
  bool
  operator==(const LeafKey& o) const {
    return a == o.a && b == o.b && c == o.c;
  }
};
struct LeafKeyHash {
  size_t
  operator()(const LeafKey& k) const {
    size_t h = std::hash<int>()(k.a);
    h = h * 1000003u ^ std::hash<int>()(k.b);
    h = h * 1000003u ^ std::hash<int>()(k.c);
    return h;
  }
};

struct Bucket {
  std::vector<std::pair<int, unsigned>> xor_cands;// (node var, raw truth)
  std::vector<std::pair<int, unsigned>> maj_cands;
};

// Enumerates all k<=3-leaf cuts of every AND node and buckets nodes whose
// cut function matches a plain XOR(k) or AND2/MAJ3(k) pattern by leaf set.
void
generate_candidates(std::vector<AdderCandidate>& out) {
  ensure_aig_view();
  const unsigned M = g_view.M;

  std::vector<std::vector<Cut>> cuts(M);
  for(unsigned v = 0; v < M; v++) {
    if(g_view.is_and[v])
      continue;
    Cut self;
    self.n = 1;
    self.leaves[0] = static_cast<int>(v);
    self.truth = base_pattern(0, 1);
    cuts[v].push_back(self);
  }

  int zero_phase[3] = { 0, 0, 0 };
  const unsigned xor2_canon = xor_pattern(2, zero_phase);
  const unsigned xor3_canon = xor_pattern(3, zero_phase);
  const unsigned mask2 = 0xFu, mask3 = 0xFFu;
  const auto& maj2 = maj_table(2);
  const auto& maj3 = maj_table(3);

  std::unordered_map<LeafKey, Bucket, LeafKeyHash> buckets;

  for(unsigned v = 0; v < M; v++) {
    if(!g_view.is_and[v])
      continue;
    const int f0 = g_view.fanin0[v], f1 = g_view.fanin1[v];
    const bool c0 = g_view.compl0[v], c1 = g_view.compl1[v];
    std::vector<Cut>& mine = cuts[v];

    for(const Cut& cut0 : cuts[f0]) {
      for(const Cut& cut1 : cuts[f1]) {
        Cut merged;
        if(!merge_cuts(cut0, c0, cut1, c1, merged))
          continue;

        bool dup = false;
        for(const Cut& existing : mine)
          if(cut_equal(existing, merged)) {
            dup = true;
            break;
          }
        if(dup)
          continue;
        mine.push_back(merged);

        if(merged.n == 2) {
          LeafKey key{ merged.leaves[0], merged.leaves[1], 0 };
          if(merged.truth == xor2_canon
             || merged.truth == ((~xor2_canon) & mask2))
            buckets[key].xor_cands.emplace_back(static_cast<int>(v),
                                                 merged.truth);
          else if(maj2.count(merged.truth))
            buckets[key].maj_cands.emplace_back(static_cast<int>(v),
                                                 merged.truth);
        } else if(merged.n == 3) {
          LeafKey key{ merged.leaves[0], merged.leaves[1], merged.leaves[2] };
          if(merged.truth == xor3_canon
             || merged.truth == ((~xor3_canon) & mask3))
            buckets[key].xor_cands.emplace_back(static_cast<int>(v),
                                                 merged.truth);
          else if(maj3.count(merged.truth))
            buckets[key].maj_cands.emplace_back(static_cast<int>(v),
                                                 merged.truth);
        }
      }
    }

    // Trivial self-cut, so this node can also serve as an opaque leaf for
    // its own consumers.
    Cut self;
    self.n = 1;
    self.leaves[0] = static_cast<int>(v);
    self.truth = base_pattern(0, 1);
    mine.push_back(self);
  }

  for(auto& kv : buckets) {
    const LeafKey& key = kv.first;
    Bucket& bucket = kv.second;
    if(bucket.xor_cands.empty() || bucket.maj_cands.empty())
      continue;

    const int n_in = (key.c != 0) ? 3 : 2;
    const int leaves[3] = { key.a, key.b, key.c };
    const auto& table = (n_in == 2) ? maj2 : maj3;
    const unsigned mask = (n_in == 2) ? mask2 : mask3;

    for(const auto& xc : bucket.xor_cands) {
      for(const auto& mc : bucket.maj_cands) {
        auto it = table.find(mc.second);
        if(it == table.end())
          continue;// unreachable: mc was only bucketed because it matched
        const PhaseInfo& mp = it->second;

        const unsigned expected_xor = xor_pattern(n_in, mp.input_phase);
        int out_phase_xor;
        if(xc.second == expected_xor)
          out_phase_xor = 0;
        else if(xc.second == ((~expected_xor) & mask))
          out_phase_xor = 1;
        else
          continue;// inconsistent pairing: not a genuine adder, skip

        AdderCandidate cand;
        cand.n_in = n_in;
        for(int i = 0; i < n_in; i++) {
          cand.in[i] = leaves[i];
          cand.input_phase[i] = mp.input_phase[i];
        }
        cand.xor_root = xc.first;
        cand.maj_root = mc.first;
        cand.output_phase_xor = out_phase_xor;
        cand.output_phase_maj = mp.output_phase;
        cand.xor_truth = xc.second;
        cand.maj_truth = mc.second;
        out.push_back(cand);
      }
    }
  }
}

// Drops (a) half-adders whose AND gate is merely a fanin of their own XOR
// gate with no other fanout (not a genuine standalone half-adder), and
// (b) half-adders fully contained inside an already-detected full-adder's
// structure (avoids double-counting a piece of an FA as its own HA).
void
filter_candidates(std::vector<AdderCandidate>& cands) {
  std::vector<AdderCandidate> kept;
  kept.reserve(cands.size());
  for(const auto& c : cands) {
    if(c.n_in == 2) {
      const bool maj_is_fanin_of_xor
        = (g_view.fanin0[c.xor_root] == c.maj_root)
          || (g_view.fanin1[c.xor_root] == c.maj_root);
      if(maj_is_fanin_of_xor && g_view.refcount[c.maj_root] == 1)
        continue;
    }
    kept.push_back(c);
  }

  std::vector<uint8_t> insider(g_view.M, 0);
  std::vector<uint8_t> visited(g_view.M, 0);
  for(const auto& c : kept) {
    if(c.n_in != 3)
      continue;
    for(int i = 0; i < 3; i++)
      visited[c.in[i]] = 1;// stop the cone walk at the FA's own inputs
    for(int root : { c.xor_root, c.maj_root }) {
      std::vector<int> stack{ root };
      while(!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        if(visited[u])
          continue;
        visited[u] = 1;
        insider[u] = 1;
        if(g_view.is_and[u]) {
          stack.push_back(g_view.fanin0[u]);
          stack.push_back(g_view.fanin1[u]);
        }
      }
    }
  }

  std::vector<AdderCandidate> final_list;
  final_list.reserve(kept.size());
  for(const auto& c : kept) {
    if(c.n_in == 2 && insider[c.xor_root] && insider[c.maj_root])
      continue;
    final_list.push_back(c);
  }
  cands.swap(final_list);
}

/*------------------------------------------------------------------------*/
// Body-node collection (which internal AND nodes belong to a candidate's
// XOR/MAJ realization), mirroring ABC's recursiveXor DFS.
/*------------------------------------------------------------------------*/
void
collect_cone(int root,
             const int* leaves,
             int n_leaves,
             std::unordered_set<int>& visited,
             std::unordered_set<int>& body) {
  auto is_leaf = [&](int x) {
    for(int i = 0; i < n_leaves; i++)
      if(leaves[i] == x)
        return true;
    return false;
  };
  std::vector<int> stack{ root };
  while(!stack.empty()) {
    int u = stack.back();
    stack.pop_back();
    if(visited.count(u))
      continue;
    visited.insert(u);
    if(!g_view.is_and[u]) {
      body.insert(u);
      continue;
    }
    const int f0 = g_view.fanin0[u], f1 = g_view.fanin1[u];
    body.insert(f0);
    body.insert(f1);
    if(!is_leaf(f0) && !visited.count(f0))
      stack.push_back(f0);
    if(!is_leaf(f1) && !visited.count(f1))
      stack.push_back(f1);
  }
}

std::vector<int>
build_node_list(const AdderCandidate& c) {
  std::vector<int> nodes;
  nodes.reserve(8);
  for(int i = 0; i < c.n_in; i++)
    nodes.push_back(c.in[i]);

  std::unordered_set<int> visited, body;
  collect_cone(c.xor_root, c.in, c.n_in, visited, body);
  collect_cone(c.maj_root, c.in, c.n_in, visited, body);

  std::vector<int> mid(body.begin(), body.end());
  std::sort(mid.begin(), mid.end());
  for(int id : mid) {
    if(id == c.xor_root || id == c.maj_root)
      continue;
    bool is_input = false;
    for(int i = 0; i < c.n_in; i++)
      if(c.in[i] == id) {
        is_input = true;
        break;
      }
    if(is_input)
      continue;
    nodes.push_back(id);
  }
  nodes.push_back(c.xor_root);
  nodes.push_back(c.maj_root);
  return nodes;
}

// The XOR/MAJ phase pairing above accepts any XOR(k) root paired with any
// AND2/MAJ3 root over the same leaves, because 2*maj + xor = a + b (+ c)
// holds for *any* such pair -- the relation is valid even when the two
// roots belong to two different physical adder cells that happen to share
// their inputs. Those coincidental pairs verify fine, but substituting them
// merges the cone ids of unrelated cells (mark_one_faha_instance()), which
// scrambles the term order the reduction relies on.
//
// A real full-adder cell realizes its carry more cheaply than its sum: the
// MAJ3 cone is strictly smaller than the XOR3 cone. A cross-paired MAJ3
// comes from somewhere else and never is. Checked on the aoki 64-bit
// benchmarks, this reproduces ABC's accept/reject decision for every one of
// the 16540 full-adder candidates. Half-adders are left to the dedupe pass
// below; for them the same measure is necessary but not sufficient.
size_t
cone_size(int root, const int* leaves, int n_leaves) {
  std::unordered_set<int> visited, body;
  collect_cone(root, leaves, n_leaves, visited, body);
  return body.size();
}

void
drop_cross_paired_full_adders(std::vector<AdderCandidate>& cands) {
  std::vector<AdderCandidate> kept;
  kept.reserve(cands.size());
  for(const auto& c : cands) {
    if(c.n_in == 3
       && cone_size(c.maj_root, c.in, c.n_in)
            >= cone_size(c.xor_root, c.in, c.n_in))
      continue;
    kept.push_back(c);
  }
  cands.swap(kept);
}

/*------------------------------------------------------------------------*/
// A single physical gate can be a locally-consistent XOR/AND(MAJ) root or
// body node for more than one candidate at once: cuts are enumerated
// exhaustively, so a node buried inside one adder's realization can also be
// the trivial leaf/root of an unrelated, coincidentally-matching candidate
// elsewhere. Each is individually a true relation, but claiming the same
// gates for many overlapping "adders" starves the reverse K-cover of a
// useful boundary and floods the linear-polynomial cache with redundant
// entries. Keep candidates non-overlapping in their internal footprint
// (everything but their leaf inputs, which legitimately fan out to many
// consumers): process full-adders before half-adders, smaller/simpler
// realizations first, and skip any candidate whose footprint has already
// been claimed.
void
dedupe_overlapping_candidates(std::vector<AdderCandidate>& cands,
                               std::vector<std::vector<int>>& node_lists) {
  std::vector<std::vector<int>> lists(cands.size());
  for(size_t i = 0; i < cands.size(); i++)
    lists[i] = build_node_list(cands[i]);

  std::vector<size_t> order(cands.size());
  for(size_t i = 0; i < order.size(); i++)
    order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if(cands[a].n_in != cands[b].n_in)
      return cands[a].n_in > cands[b].n_in;// full-adders before half-adders
    if(lists[a].size() != lists[b].size())
      return lists[a].size() < lists[b].size();// simpler realizations first
    return cands[a].xor_root < cands[b].xor_root;// deterministic tie-break
  });

  auto is_input = [](const AdderCandidate& c, int id) {
    for(int i = 0; i < c.n_in; i++)
      if(c.in[i] == id)
        return true;
    return false;
  };

  std::vector<uint8_t> claimed(g_view.M, 0);
  std::vector<AdderCandidate> kept;
  std::vector<std::vector<int>> kept_lists;
  kept.reserve(cands.size());
  kept_lists.reserve(cands.size());
  for(size_t idx : order) {
    const AdderCandidate& c = cands[idx];
    const std::vector<int>& nodes = lists[idx];

    bool overlap = false;
    for(int id : nodes)
      if(!is_input(c, id) && claimed[id]) {
        overlap = true;
        break;
      }
    if(overlap)
      continue;

    for(int id : nodes)
      if(!is_input(c, id))
        claimed[id] = 1;
    kept.push_back(c);
    kept_lists.push_back(nodes);
  }
  cands.swap(kept);
  node_lists.swap(kept_lists);
}

void
populate_structs(const std::vector<AdderCandidate>& cands,
                  const std::vector<std::vector<int>>& node_lists,
                  std::vector<HalfAdder>& has,
                  std::vector<FullAdder>& fas) {
  int next_ha_id = 0, next_fa_id = 0;
  for(size_t idx = 0; idx < cands.size(); idx++) {
    const AdderCandidate& c = cands[idx];
    std::vector<int> nodes = node_lists[idx];
    std::vector<uint8_t> tt
      = { static_cast<uint8_t>(c.xor_truth), static_cast<uint8_t>(c.maj_truth) };
    std::vector<int> out_phase = { c.output_phase_xor, c.output_phase_maj };

    if(c.n_in == 2) {
      HalfAdder h;
      h.id = next_ha_id++;
      h.nodes = std::move(nodes);
      h.tt = std::move(tt);
      h.input_phase = { c.input_phase[0], c.input_phase[1] };
      h.output_phase = std::move(out_phase);
      has.push_back(std::move(h));
    } else {
      FullAdder f;
      f.id = next_fa_id++;
      f.nodes = std::move(nodes);
      f.tt = std::move(tt);
      f.input_phase
        = { c.input_phase[0], c.input_phase[1], c.input_phase[2] };
      f.output_phase = std::move(out_phase);
      fas.push_back(std::move(f));
    }
  }
}

/*------------------------------------------------------------------------*/
// dedupe_overlapping_candidates() exists to stop two *different* adders from
// claiming the same gate. A duplicated realization of the *same* adder --
// same leaves, sharing one of its two roots, differing only in which of
// several equivalent gates computes the other output -- is not a competing
// adder, and the K-cover needs to see all of them: which root a primary
// output meets first decides the boundary it gets. Put those back, but
// nothing else.
std::vector<AdderCandidate>
augment_with_duplicate_realizations(const std::vector<AdderCandidate>& kept,
                                    const std::vector<AdderCandidate>& pool) {
  auto leaf_key = [](const AdderCandidate& c) {
    return LeafKey{ c.in[0], c.in[1], c.n_in == 3 ? c.in[2] : 0 };
  };
  std::unordered_map<LeafKey, std::vector<const AdderCandidate*>, LeafKeyHash>
    by_leaves;
  for(const auto& c : kept)
    by_leaves[leaf_key(c)].push_back(&c);

  std::vector<AdderCandidate> out = kept;
  for(const auto& c : pool) {
    auto it = by_leaves.find(leaf_key(c));
    if(it == by_leaves.end()) {
      // dedupe dropped every realization of this cut, because some other
      // adder had already claimed the gates. ABC's pool never drops a
      // detected cut, and the boundary search needs to see it, so put it
      // back.
      if(c.n_in == 3)
        out.push_back(c);
      continue;
    }
    for(const AdderCandidate* k : it->second) {
      if(k->xor_root == c.xor_root && k->maj_root == c.maj_root)
        break;// already kept
      if(k->xor_root == c.xor_root || k->maj_root == c.maj_root) {
        out.push_back(c);
        break;
      }
    }
  }
  return out;
}

/*------------------------------------------------------------------------*/
// Cache of the adder-pair list, so get_kcover_nodes() can reuse what
// detect_ha_fa_xor() already computed.
//
// The K-cover gets its own view of the candidates: see
// augment_with_duplicate_realizations(). Handing it the substitution list
// makes some primary outputs latch onto a deeper adder pair, leaving the
// final-stage adder underconstrained; handing it the completely unthinned
// pool goes too far the other way and cuts the boundary short.
/*------------------------------------------------------------------------*/
struct DetectionCache {
  bool valid = false;
  std::vector<AdderCandidate> kcover_pool;
};
DetectionCache g_cache;

/*------------------------------------------------------------------------*/
// Reverse local K-cover: for each primary output, find the shallowest
// reachable adder pair and grow a local cone bounded by that pair's inputs
// and by AIG level.
/*------------------------------------------------------------------------*/
struct PairMaps {
  std::vector<int> root_to_pair;// size M, -1 if var is not an adder root
  std::vector<std::vector<int>> pair_inputs;
  std::vector<std::vector<int>> pair_roots;
};

PairMaps
build_pair_maps(const std::vector<AdderCandidate>& adders) {
  PairMaps m;
  m.root_to_pair.assign(g_view.M, -1);
  m.pair_inputs.resize(adders.size());
  m.pair_roots.resize(adders.size());
  for(size_t i = 0; i < adders.size(); i++) {
    const auto& c = adders[i];
    for(int k = 0; k < c.n_in; k++)
      m.pair_inputs[i].push_back(c.in[k]);
    m.pair_roots[i] = { c.xor_root, c.maj_root };
    m.root_to_pair[c.xor_root] = static_cast<int>(i);
    m.root_to_pair[c.maj_root] = static_cast<int>(i);
  }
  return m;
}

// Level-order BFS from start_var, returning the shallowest un-banned adder
// root reachable; prefers a pair not yet claimed by an earlier PO over one
// that already was.
void
find_first_hit_pair(int start_var,
                     const PairMaps& maps,
                     const std::vector<uint8_t>& pair_hit,
                     const std::vector<uint8_t>& banned,
                     int& out_root,
                     int& out_pair) {
  out_root = -1;
  out_pair = -1;
  if(start_var <= 0)
    return;

  std::vector<int> q{ start_var };
  std::vector<uint8_t> visited(g_view.M, 0);
  visited[start_var] = 1;

  size_t head = 0;
  while(head < q.size()) {
    const size_t level_begin = head, level_end = q.size();
    int best_new_root = -1, best_new_pair = -1;
    int best_old_root = -1, best_old_pair = -1;

    for(size_t idx = level_begin; idx < level_end; idx++) {
      const int u = q[idx];
      if(banned[u])
        continue;
      const int pid = maps.root_to_pair[u];
      if(pid < 0)
        continue;
      if(!pair_hit[pid]) {
        if(best_new_pair < 0) {
          best_new_root = u;
          best_new_pair = pid;
        }
      } else if(best_old_pair < 0) {
        best_old_root = u;
        best_old_pair = pid;
      }
    }
    if(best_new_pair >= 0) {
      out_root = best_new_root;
      out_pair = best_new_pair;
      return;
    }
    if(best_old_pair >= 0) {
      out_root = best_old_root;
      out_pair = best_old_pair;
      return;
    }

    for(size_t idx = level_begin; idx < level_end; idx++) {
      const int u = q[idx];
      if(!g_view.is_and[u])
        continue;
      // Visit fanins in increasing node id. Which of two roots sitting on
      // the same BFS level is reported first decides the boundary this PO
      // gets, so the order has to be pinned down. AIGER stores an AND as
      // rhs0 >= rhs1, i.e. fanin0 is the *larger* id, so iterate fanin1
      // before fanin0.
      const int lo = std::min(g_view.fanin0[u], g_view.fanin1[u]);
      const int hi = std::max(g_view.fanin0[u], g_view.fanin1[u]);
      for(int f : { lo, hi })
        if(f > 0 && !visited[f]) {
          visited[f] = 1;
          q.push_back(f);
        }
    }
    head = level_end;
  }
}

// DFS from 'roots' through AND fanins, stopping at nodes already marked as
// boundary inputs or at PI/const; marks every visited node (including stop
// points) into 'hit'.
void
collect_pair_cone(const std::vector<int>& roots,
                   const std::vector<uint8_t>& boundary_inputs,
                   std::vector<uint8_t>& hit) {
  std::vector<int> stack;
  for(int r : roots)
    if(r > 0)
      stack.push_back(r);
  while(!stack.empty()) {
    int u = stack.back();
    stack.pop_back();
    if(u <= 0 || hit[u])
      continue;
    hit[u] = 1;
    if(boundary_inputs[u] || !g_view.is_and[u])
      continue;
    stack.push_back(g_view.fanin0[u]);
    stack.push_back(g_view.fanin1[u]);
  }
}

}// namespace

/*------------------------------------------------------------------------*/
bool
detect_ha_fa_xor(std::vector<HalfAdder>& half_adders,
                  std::vector<FullAdder>& full_adders,
                  std::vector<XorRemain>& xor_list) {
  half_adders.clear();
  full_adders.clear();
  xor_list.clear();// leftover XOR relations are already covered by affine
                    // cuts elsewhere; kept for API compatibility only

  aiger* model = get_model();
  if(!model) {
    std::cerr << "Error: AIG has not been parsed" << std::endl;
    return false;
  }
  ensure_aig_view();

  std::vector<AdderCandidate> raw;
  generate_candidates(raw);
  filter_candidates(raw);

  const std::vector<AdderCandidate> pool = raw;// before either pass below

  // The K-cover and the substitution list want different things out of the
  // same candidates, so they are derived separately.
  //
  // K-cover: dedupe, then put back the entries dedupe cost us. It must not
  // see the cross-paired full adders removed either -- ABC likewise runs its
  // boundary search over the raw detected-adder pool.
  std::vector<AdderCandidate> kc_base = pool;
  {
    std::vector<std::vector<int>> kc_lists;
    dedupe_overlapping_candidates(kc_base, kc_lists);
  }
  g_cache.kcover_pool = augment_with_duplicate_realizations(kc_base, pool);
  g_cache.valid = true;

  // Substitution list: drop the coincidental pairings, then dedupe so no
  // gate is claimed by two adders.
  drop_cross_paired_full_adders(raw);

  std::vector<std::vector<int>> node_lists;
  dedupe_overlapping_candidates(raw, node_lists);

  populate_structs(raw, node_lists, half_adders, full_adders);

  if(verbose >= 1)
    std::cout << "Detected: half_adders=" << half_adders.size()
              << ", full_adders=" << full_adders.size()
              << ", xor_remaining=" << xor_list.size() << std::endl;

  return true;
}

/*------------------------------------------------------------------------*/
bool
get_kcover_nodes(std::vector<int>& reverse_fsa_nodes,
                 std::vector<int>& boundary_nodes) {
  reverse_fsa_nodes.clear();
  boundary_nodes.clear();

  aiger* model = get_model();
  if(!model) {
    std::cerr << "Error: AIG has not been parsed" << std::endl;
    return false;
  }
  ensure_aig_view();

  std::vector<AdderCandidate> local_adders;
  const std::vector<AdderCandidate>* adders_ptr;
  if(g_cache.valid) {
    if(verbose >= 1)
      std::cout << "Using cached adder pairs for K-Cover computation"
                << std::endl;
    adders_ptr = &g_cache.kcover_pool;
  } else {
    generate_candidates(local_adders);
    filter_candidates(local_adders);
    adders_ptr = &local_adders;
  }
  const std::vector<AdderCandidate>& adders = *adders_ptr;

  PairMaps maps = build_pair_maps(adders);
  std::vector<uint8_t> pair_hit(adders.size(), 0);
  std::vector<uint8_t> banned(g_view.M, 0);
  std::vector<uint8_t> boundary_inputs(g_view.M, 0);

  const unsigned num_po = model->num_outputs;
  std::vector<std::vector<int>> po_inputs(num_po);
  bool any_boundary = false;

  for(unsigned i = 0; i < num_po; i++) {
    const int var = static_cast<int>(model->outputs[i].lit >> 1);
    if(var <= 0)
      continue;

    int root_id = -1, pair_id = -1;
    find_first_hit_pair(var, maps, pair_hit, banned, root_id, pair_id);
    if(pair_id < 0)
      continue;

    pair_hit[pair_id] = 1;
    for(int x : maps.pair_inputs[pair_id])
      if(x > 0) {
        boundary_inputs[x] = 1;
        any_boundary = true;
      }

    std::vector<int> slot = maps.pair_inputs[pair_id];
    std::sort(slot.begin(), slot.end());
    slot.erase(std::unique(slot.begin(), slot.end()), slot.end());
    po_inputs[i] = std::move(slot);

    std::vector<uint8_t> tmp_hit(g_view.M, 0);
    collect_pair_cone(maps.pair_roots[pair_id], boundary_inputs, tmp_hit);
    for(unsigned j = 0; j < g_view.M; j++)
      if(tmp_hit[j])
        banned[j] = 1;
  }

  // Grow, per PO (in output order, so later POs see earlier POs' borders),
  // a local cone from the PO down to its recorded adder inputs / the
  // circuit's global border / PIs, additionally cut off by AIG level.
  std::vector<uint8_t> border_all(g_view.M, 0);
  std::vector<uint8_t> middle_all(g_view.M, 0);

  for(unsigned i = 0; i < num_po; i++) {
    const int var = static_cast<int>(model->outputs[i].lit >> 1);
    if(var <= 0)
      continue;

    const std::vector<int>& po_border_vec = po_inputs[i];
    std::vector<uint8_t> po_border(g_view.M, 0);
    for(int x : po_border_vec)
      po_border[x] = 1;

    const bool has_min_lv = !po_border_vec.empty();
    int min_lv = 0;
    if(has_min_lv) {
      min_lv = g_view.level[po_border_vec[0]];
      for(int x : po_border_vec)
        min_lv = std::min(min_lv, g_view.level[x]);
    }

    std::vector<int> local_middle, new_border;
    std::vector<uint8_t> visited(g_view.M, 0);
    std::vector<int> stack{ var };
    while(!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      if(u <= 0 || visited[u])
        continue;
      visited[u] = 1;

      const bool is_po_border = po_border[u];
      const bool is_glob_border = border_all[u];
      const bool is_pi_or_const = !g_view.is_and[u];

      if(is_po_border || is_glob_border || is_pi_or_const) {
        if(!is_glob_border)
          new_border.push_back(u);
        continue;
      }

      local_middle.push_back(u);
      if(has_min_lv && g_view.level[u] < min_lv) {
        new_border.push_back(u);
        continue;
      }
      stack.push_back(g_view.fanin0[u]);
      stack.push_back(g_view.fanin1[u]);
    }

    for(int u : new_border)
      border_all[u] = 1;
    for(int u : local_middle)
      if(!border_all[u])
        middle_all[u] = 1;
    for(int u : new_border)
      middle_all[u] = 0;// border takes priority on conflict
  }

  // Closure: every fanin of a middle node must itself be middle or border.
  bool changed = true;
  while(changed) {
    changed = false;
    std::vector<int> round_border;
    for(unsigned id = 1; id < g_view.M; id++) {
      if(!middle_all[id] || !g_view.is_and[id])
        continue;
      for(int f : { g_view.fanin0[id], g_view.fanin1[id] })
        if(!middle_all[f] && !border_all[f]) {
          border_all[f] = 1;
          round_border.push_back(f);
          changed = true;
        }
    }
    for(int f : round_border)
      middle_all[f] = 0;
  }

  for(unsigned id = 1; id < g_view.M; id++)
    if(border_all[id] || middle_all[id])
      reverse_fsa_nodes.push_back(static_cast<int>(id));
  for(unsigned id = 1; id < g_view.M; id++)
    if(border_all[id])
      boundary_nodes.push_back(static_cast<int>(id));

  // No adder pair fed any PO: no useful boundary could be established, so
  // treat the whole AIG as free simplification area.
  if(!any_boundary) {
    reverse_fsa_nodes.clear();
    for(unsigned id = 0; id < g_view.M; id++)
      reverse_fsa_nodes.push_back(static_cast<int>(id));
  }

  if(verbose >= 1)
    std::cout << "K-Cover: reverse_fsa_nodes=" << reverse_fsa_nodes.size()
              << ", boundary_nodes=" << boundary_nodes.size() << std::endl;

  return true;
}

/*------------------------------------------------------------------------*/
void
clear_detection_cache() {
  g_cache.valid = false;
  g_cache.kcover_pool.clear();
  g_cache.kcover_pool.shrink_to_fit();
}

/*------------------------------------------------------------------------*/
void
print_aig_stats() {
  if(verbose < 1)
    return;
  aiger* model = get_model();
  if(!model)
    return;
  std::cout << "=== AIG Statistics ===" << std::endl;
  std::cout << "  Primary Inputs (PIs):  " << model->num_inputs << std::endl;
  std::cout << "  Primary Outputs (POs): " << model->num_outputs
            << std::endl;
  std::cout << "  AND gates:             " << model->num_ands << std::endl;
  std::cout << "  Total objects:         "
            << (1 + model->num_inputs + model->num_ands) << std::endl;
  std::cout << "========================" << std::endl;
}

}// namespace talisman
