/*------------------------------------------------------------------------*/
/*! \file linear_algebra.cpp
    \brief contains linear algebra

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "linear_algebra.hpp"

#include <cassert>
#include <chrono>
#include <unordered_map>

#include "fast_mod.hpp"
#include "nmod.h"
#include "nmod_vec.h"

using namespace talisman;

namespace talisman {

std::vector<uint64_t>
make_bit_pattern_batch(int n, uint64_t batch_index) {
  std::vector<uint64_t> result(n, 0);

  uint64_t base = batch_index * 64;// starting row for this batch
  uint64_t rows_in_this_batch = std::min(64ULL, (1ULL << n) - base);

  for(int col = 0; col < n; ++col) {
    uint64_t word = 0;
    for(uint64_t bit = 0; bit < rows_in_this_batch; ++bit) {
      uint64_t row = base + bit;

      // Get bit 'col' of row
      if((row >> col) & 1ULL) {
        // Set the bit counting from the left (highest bit = row 0)
        word |= (1ULL << (63 - bit));
      }
    }
    result[col] = word;
  }
  return result;
}
/*------------------------------------------------------------------------*/
// assumes that all zero rows are at the bottom of the matrix
slong
first_zero_row(nmod_mat_t mat) {
  slong lo = 0;
  slong hi = mat->r;

  while(lo < hi) {
    slong mid = lo + (hi - lo) / 2;
    if(nmod_mat_is_zero_row(mat, mid))
      hi = mid;
    else
      lo = mid + 1;
  }
  return lo;
}
/*------------------------------------------------------------------------*/

void
append_repair_samples(nmod_mat_t mat, std::vector<std::vector<int>>& samples) {

  slong start = first_zero_row(mat);

  if(start + samples.size() > mat->r)
    die(44, "Running out of rows in repair sample append");

  for(size_t i = 0; i < samples.size(); i++) {
    slong row_idx = start + i;
    ulong* row = nmod_mat_entry_ptr(mat, row_idx, 0);

    const std::vector<int>& sample = samples[i];

    for(size_t j = 0; j < sample.size(); j++)
      row[j] = sample[j];
  }
}
/*------------------------------------------------------------------------*/

// compute mat[i] - a*mat[j]
void
xmay(nmod_mat_t mat, int i, int j, int a) {

  ulong* row_i = nmod_mat_entry_ptr(mat, i, 0);
  ulong* row_j = nmod_mat_entry_ptr(mat, j, 0);
  _nmod_vec_scalar_addmul_nmod(row_i, row_j, mat->c, mat->mod.n - a, mat->mod);
}
/*------------------------------------------------------------------------*/
// normalize row i starting from index j
inline void
normalize_row(ulong* row, int j, slong cols, nmod_t mod) {
  ulong c = row[j];
  if(c == 1)
    return;

  ulong inv = nmod_inv(c, mod);
  row[j] = 1;
  _nmod_vec_scalar_mul_nmod_shoup(
    row + j + 1, row + j + 1, cols - j - 1, inv, mod);
}
/*------------------------------------------------------------------------*/
slong
permute_rows(nmod_mat_t mat, std::vector<slong>& piv) {
  std::vector<slong> perm_act(mat->r);
  slong next_row = 0;
  slong rank = 0;

  // pivot rows first
  for(slong j = 0; j < mat->c; j++) {
    if(piv[j] >= 0) {
      perm_act[next_row++] = piv[j];
      rank++;
    }
  }

  // remaining rows
  for(slong i = 0; i < mat->r; i++) {
    if(std::find(piv.begin(), piv.end(), i) == piv.end()) {
      perm_act[next_row++] = i;
    }
  }

  nmod_mat_permute_rows(mat, perm_act.data(), NULL);
  std::vector<slong> perm_inv(mat->r);
  for(slong i = 0; i < mat->r; i++) {
    perm_inv[perm_act[i]] = i;
  }

  for(slong j = 0; j < mat->c; j++) {
    if(piv[j] >= 0) {
      piv[j] = perm_inv[piv[j]];
    }
  }
  return rank;
}

/*------------------------------------------------------------------------*/
void
backward_elimination(nmod_mat_t mat, std::vector<slong>& piv) {
  // backward elimination
  for(slong j = mat->c - 1; j >= 0; j--) {
    slong pivot_row = piv[j];
    if(pivot_row < 0)
      continue;

    ulong* row_p = nmod_mat_entry_ptr(mat, pivot_row, 0);

    for(slong i = 0; i < pivot_row; i++) {
      ulong* row_i = nmod_mat_entry_ptr(mat, i, 0);
      ulong val = row_i[j];
      if(val != 0)
        _nmod_vec_scalar_addmul_nmod(
          row_i + j, row_p + j, mat->c - j, mat->mod.n - val, mat->mod);
    }
  }
}
/*------------------------------------------------------------------------*/
inline bool
reduce_row(nmod_mat_t mat, ulong* row, std::vector<slong>& piv, int i) {

  const nmod_t mod = mat->mod;
  const slong cols = mat->c;

  for(slong j = 0; j < cols; j++) {
    ulong val = row[j];
    if(val == 0)
      continue;

    slong p = piv[j];
    if(p < 0) {
      normalize_row(row, j, cols, mod);
      piv[j] = i;
      return false;
    }
    ulong* row_p = nmod_mat_entry_ptr(mat, p, 0);
    _nmod_vec_scalar_addmul_nmod(
      row + j, row_p + j, cols - j, nmod_neg(val, mod), mod);
  }
  return true;
}
/*------------------------------------------------------------------------*/
slong
rref(nmod_mat_t mat, int nr_appened_samples) {

  slong relevant_rows = first_zero_row(mat);
  slong rank = relevant_rows - nr_appened_samples;

  std::vector<slong> piv(mat->c, -1);
  // collect pivots
  for(slong i = 0; i < rank; ++i) {
    ulong* row = nmod_mat_entry_ptr(mat, i, 0);
    for(slong j = 0; j < mat->c; ++j) {
      if(row[j] != 0) {
        piv[j] = i;
        break;
      }
    }
  }

  // forward elimination
  for(slong i = rank; i < relevant_rows; i++) {
    ulong* row_i = nmod_mat_entry_ptr(mat, i, 0);
    reduce_row(mat, row_i, piv, i);
  }

  // permute rows
  rank = permute_rows(mat, piv);

  // backward elimination
  backward_elimination(mat, piv);

  return rank;
}
/*------------------------------------------------------------------------*/
slong
randomized_rref(nmod_mat_t mat, slong new_rows) {
  const slong rows = mat->r;
  const slong cols = mat->c;
  const nmod_t mod = mat->mod;

  slong rank = 0;

  // copy to bit-matrix data structure
  bit_matrix<false> bitmat;
  int bitrows = 64 * (rows / 64 + 1);
  bitmat.init(bitrows, cols);
  bitmat.zero();
  for(slong i = 0; i < rows; i++) {
    ulong* row = nmod_mat_entry_ptr(mat, i, 0);
    for(slong j = 0; j < cols; j++) {
      if(row[j] != 0)
        bitmat.set_entry(i, j, 1);
    }
  }

  // make random linear combinations of rows
  for(int i = 0; i < new_rows; i++) {
    ulong* row = nmod_mat_entry_ptr(mat, i, 0);
    bitmat.random_combination(row);
    _nmod_vec_reduce(row, row, cols, mod);
  }

  // compute rref of new matrix
  nmod_mat_t upper;
  nmod_mat_window_init(upper, mat, 0, 0, new_rows, cols);
  rank = nmod_mat_rref(upper);
  nmod_mat_window_clear(upper);

  // collect pivots
  std::vector<slong> piv(mat->c, -1);
  for(slong i = 0; i < new_rows; i++) {
    ulong* row_i = nmod_mat_entry_ptr(mat, i, 0);
    for(slong j = 0; j < mat->c; j++) {
      if(row_i[j] != 0) {
        piv[j] = i;
        break;
      }
    }
  }

  // test random linear combinations
  int failures = 0;
  slong next_row = new_rows;
  bool is_zero;
  int max_rank = std::min(mat->r, mat->c);
  while(failures < 3 and rank < max_rank) {
    ulong* row = nmod_mat_entry_ptr(mat, next_row, 0);
    bitmat.random_combination(row);
    _nmod_vec_reduce(row, row, cols, mod);
    is_zero = reduce_row(mat, row, piv, next_row);
    if(is_zero)
      failures++;
    else {
      failures = 0;
      next_row++;
      rank++;
    }
  }

  // permute rows correctly, just to be safe
  // removing this has caused problems on linux
  rank = permute_rows(mat, piv);

  // if we have new pivots -> backward elimination
  if(next_row != new_rows)
    backward_elimination(mat, piv);

  // fill bottom with zeros
  for(int i = next_row; i < rows; i++) {
    ulong* row = nmod_mat_entry_ptr(mat, i, 0);
    std::fill(row, row + cols, 0);
  }

  // release the temporary bit-matrix storage before returning
  bitmat.clear();

  return rank;
}

int
nonzero_entries_in_row(ulong* row, int n) {
  int nnz = 0;
  for(int i = 0; i < n; i++)
    if(row[i] != 0)
      nnz++;
  return nnz;
}

thread_local std::unordered_map<ulong, int> counts;
bool
sparsify_row(ulong* row_i, ulong* row_j, int n, nmod_t mod) {

  int new_nnz = 0;
  ulong best_c = 0;
  int best_count = 0;

  counts.clear();
  counts.reserve(n);

  for(int k = 0; k < n; k++) {

    new_nnz += row_i[k] == 0 and row_j[k] != 0 ? 1 : 0;
    if(row_i[k] == 0 || row_j[k] == 0)
      continue;

    ulong c = nmod_div(mod.n - row_i[k], row_j[k], mod);

    int count = ++counts[c];
    if(count > best_count) {
      best_count = count;
      best_c = c;
    }
  }

  ulong p = mod.n;
  if(best_count > new_nnz)
    for(int k = 0; k < n; k++)
      row_i[k] = (row_i[k] + best_c * row_j[k]) % p;

  return best_count > new_nnz;
}

bool
sparsify(nmod_mat_t mat) {
  int c = mat->c;
  bool change = false;
  for(int i = mat->r - 1; i >= 0; i--) {
    auto row_i = nmod_mat_entry_ptr(mat, i, 0);
    if(nonzero_entries_in_row(row_i, c) < 8)
      continue;
    for(int j = i - 1; j >= 0; j--) {
      change |= sparsify_row(row_i, nmod_mat_entry_ptr(mat, j, 0), c, mat->mod);
    }
  }
  return change;
}
/*------------------------------------------------------------------------*/
thread_local std::vector<slong> pivots;
thread_local std::vector<slong> nonpivots;
extern bool custom_linalg;
size_t
kernel(nmod_mat_t mat, nmod_mat_t K, int nr_appended_samples) {
  slong i, j, k, n, nullity, rank;

  n = mat->c;

  if(nr_appended_samples > 0 and custom_linalg)
    rank = rref(mat, nr_appended_samples);
  else if(mat->c > 100 and custom_linalg)
    rank = randomized_rref(mat, 2 * mat->c / 3);
  else
    rank = nmod_mat_rref(mat);

  nullity = n - rank;

  pivots.clear();
  nonpivots.clear();
  pivots.reserve(rank);
  nonpivots.reserve(nullity);

  nmod_mat_init(K, nullity, n, mat->mod.n);
  nmod_mat_zero(K);

  if(rank == 0) {
    for(i = 0; i < nullity; i++)
      nmod_mat_entry(K, i, i) = 1;
  } else if(nullity) {
    for(i = j = k = 0; i < rank; i++) {
      while(nmod_mat_entry(mat, i, j) == 0) {
        nonpivots[k] = j;
        k++;
        j++;
      }
      pivots[i] = j;
      j++;
    }
    while(k < nullity) {
      nonpivots[k] = j;
      k++;
      j++;
    }

    for(i = 0; i < nullity; i++) {
      for(j = 0; j < rank; j++) {
        unsigned long c = nmod_mat_entry(mat, j, nonpivots[i]);
        nmod_mat_entry(K, i, pivots[j]) = nmod_neg(c, mat->mod);
      }

      nmod_mat_entry(K, i, nonpivots[i]) = 1;
    }
  }

  if(nullity and nmod_mat_entry(K, K->r - 1, K->c - 1) != 0) {
    for(int i = 0; i < 4; i++)
      if(!sparsify(K))
        break;
  }

  return nullity;
}

/*------------------------------------------------------------------------*/

}// namespace talisman
