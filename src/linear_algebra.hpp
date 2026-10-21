/*------------------------------------------------------------------------*/
/*! \file linear_algebra.hpp
    \brief contains the linear algebra

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_LINEAR_ALGEBRA_H_
#define TALISMAN_SRC_LINEAR_ALGEBRA_H_
/*------------------------------------------------------------------------*/

#include "gate.hpp"
#include "libpopcnt.h"

#include <nmod_mat.h>
#include <random>
#include <vector>
/*------------------------------------------------------------------------*/

namespace talisman {

#define GETBIT(n, k) (((n) >> (k)) & 1)
#define SETBIT(n, k, b) \
  ((n) = ((n) & ~(1LL << (k))) | ((sample_t)(!!(b)) << (k)))

extern const size_t batch_size;

void inline set_input_matrix(nmod_mat_t mat, size_t n) {
  assert(nmod_mat_nrows(mat) >= n);
  for(size_t i = 0; i < n; i++)
    nmod_mat_set_entry(mat, i, n - i + 1, 1);
}

void inline set_constant_column(nmod_mat_t mat, size_t j) {
  for(signed long i = 0; i < nmod_mat_nrows(mat); ++i) {
    nmod_mat_set_entry(mat, i, j, 1);
  }
}

// Find indices of columns in an nmod_mat_t that are constant (all rows
// contain the same value). Returns an empty vector if no constant
// columns found.
inline std::vector<int>
find_constant_columns(nmod_mat_t M) {
  int r = nmod_mat_nrows(M);
  int c = nmod_mat_ncols(M);
  std::vector<int> out;
  if(r == 0 || c == 0)
    return out;

  for(int j = 0; j < c; ++j) {
    int base = nmod_mat_entry(M, 0, j);
    bool constant = true;
    for(int i = 1; i < (int)r; ++i) {
      if(nmod_mat_entry(M, i, j) != base) {
        constant = false;
        break;
      }
    }
    if(constant)
      out.push_back(j);
  }
  return out;
}

// insert the sample s into column j and rows from start to max(start +
// batch_size, max)
void inline mat_insert_batch(nmod_mat_t mat,
                             size_t start,
                             size_t j,
                             sample_t s,
                             size_t max) {

  for(size_t i = 0; i < batch_size && start + i < max; ++i) {
    uint64_t bit = (s >> (batch_size - 1 - i)) & 1;
    nmod_mat_set_entry(mat, start + i, j, bit);
  }
}
/*------------------------------------------------------------------------*/
std::vector<uint64_t>
make_bit_pattern_batch(int n, uint64_t batch_index);
/*------------------------------------------------------------------------*/
void
xmay(nmod_mat_t mat, int i, int j, int a);
/*------------------------------------------------------------------------*/
void
append_repair_samples(nmod_mat_t mat, std::vector<std::vector<int>>& samples);

/*------------------------------------------------------------------------*/
size_t
kernel(nmod_mat_t mat, nmod_mat_t K, int nr_appended_samples = 0);

/*------------------------------------------------------------------------*/
template<bool row_format>
class bit_matrix {
  public:
  int rows, cols, stride;
  sample_t* entries;
  std::vector<uint64_t> mask;

  bit_matrix() {};

  void init(int r, int c) {
    cols = c;
    rows = r;
    if(row_format) {
      if(cols % 64 != 0)
        die(64, "Expect divisibility by 64");
      stride = cols / 64;
      entries = new sample_t[rows * stride]();
    } else {
      if(rows % 64 != 0)
        die(64, "Expect divisibility by 64");
      stride = rows / 64;
      entries = new sample_t[cols * stride]();
    }
  }

  void clear() { delete[] entries; }

  void print() {
    std::cout << "Matrix of size " << rows << " x " << cols << std::endl;

    for(int i = 0; i < rows; ++i) {
      for(int j = 0; j < cols; ++j)
        std::cout << entry(i, j);
      std::cout << std::endl;
    }
  }

  void inline zero() {
    if(row_format)
      std::fill(entries, entries + rows * stride, 0);
    else
      std::fill(entries, entries + cols * stride, 0);
  }

  int inline entry(int i, int j) {
    if(row_format) {
      sample_t s = entries[i * stride + (j / 64)];
      return GETBIT(s, 63 - (j % 64));
    } else {
      sample_t s = entries[j * stride + (i / 64)];
      return GETBIT(s, 63 - (i % 64));
    }
  }

  sample_t& get_row_block(int i, int b) { return entries[i * stride + b]; }
  sample_t& get_col_block(int j, int b) { return entries[j * stride + b]; }

  void inline set_constant_column(int j) {
    for(int k = 0; k < stride; k++)
      set_batch_col(64 * k, j, sample_t(-1));
  }

  void inline set_entry(int i, int j, sample_t value) {
    if(row_format) {
      sample_t& s = entries[i * stride + (j / 64)];
      s = (value << (63 - (j % 64))) | s;
    } else {
      sample_t& s = entries[j * stride + (i / 64)];
      s = (value << (63 - (i % 64))) | s;
    }
  }

  void inline set_batch_col(int i, int j, sample_t s) {
    if(i % 64 != 0)
      die(4, "Expect divisbility by 64");
    entries[j * stride + (i / 64)] = s;
  }

  void inline gram_matrix(nmod_mat_t out, int32_t prime) {
    sample_t* inner = new sample_t[stride];

    for(int i = 0; i < cols; ++i) {
      sample_t* ci = entries + i * stride;
      for(int j = i; j < cols; ++j) {
        sample_t* cj = entries + j * stride;
        for(int k = 0; k < stride; ++k)
          inner[k] = ci[k] & cj[k];
        uint64_t v = popcnt(inner, sizeof(sample_t) * stride) % prime;
        nmod_mat_set_entry(out, i, j, v);
        if(i != j)
          nmod_mat_set_entry(out, j, i, v);
      }
    }

    delete[] inner;
  }

  std::mt19937_64 gen;
  std::uniform_int_distribution<uint64_t> dis;

  void random_combination(ulong* row) {
    mask.clear();
    mask.resize(stride);
    for(slong s = 0; s < stride; s++)
      mask[s] = dis(gen);

    for(slong j = 0; j < cols; j++) {
      for(slong s = 0; s < stride; s++) {
        uint64_t x = entries[j * stride + s] & mask[s];
        row[j] += __builtin_popcountll(x);
      }
    }
  }

  // Expand the number of rows by `extra_rows`. Only supported for
  // `row_format == false` (column-major packed by rows). New row count
  // is rounded up to a multiple of 64.
  void expand_rows(int extra_rows) {
    if(extra_rows <= 0)
      return;
    if(row_format)
      die(4, "expand_rows not supported for row_format=true");
    int old_rows = rows;
    int new_rows = old_rows + extra_rows;
    if(new_rows % 64 != 0)
      new_rows = ((new_rows + 63) / 64) * 64;
    int old_stride = stride;
    int new_stride = new_rows / 64;
    sample_t* new_entries = new sample_t[cols * new_stride]();
    for(int j = 0; j < cols; ++j) {
      for(int b = 0; b < old_stride; ++b)
        new_entries[j * new_stride + b] = entries[j * old_stride + b];
    }
    delete[] entries;
    entries = new_entries;
    rows = new_rows;
    stride = new_stride;
  }
};

// Find constant columns in a packed bit_matrix<false> (columns are stored
// as stride 64-bit blocks). Returns indices of columns that are identical
// across all rows.
inline std::vector<int>
find_constant_columns(bit_matrix<false>& A) {
  std::vector<int> out;
  if(A.rows == 0 || A.cols == 0)
    return out;
  int blocks = A.stride;
  for(int j = 0; j < A.cols; ++j) {
    sample_t first = A.get_col_block(j, 0);
    bool constant = true;
    for(int b = 1; b < blocks; ++b) {
      if(A.get_col_block(j, b) != first) {
        constant = false;
        break;
      }
    }
    if(constant)
      out.push_back(j);
  }
  return out;
}

}// namespace talisman

#endif// TALISMAN_SRC_LINEAR_ALGEBRA_H_
