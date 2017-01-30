/**
 * cartan_equiv.h
 * Copyright 2016 John Lawson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#ifndef REFL_CARTAN_EQUIV_H__
#define REFL_CARTAN_EQUIV_H__

#include <armadillo>
#include <boost/dynamic_bitset.hpp>

#include <deque>

namespace refl {
/**
 * Determine whether two cartan matrices are equivalent up to simultaneous
 * changes of signs in rows and columns.
 *
 * A Cartan matrix is symmetric, and this method assumes that. Only the lower
 * triangle of the matrices is checked. The upper triangle is ignored.
 */
class CartanEquiv {
 public:
  template <class elem_t>
  bool operator()(arma::Mat<elem_t> const& lhs, arma::Mat<elem_t> const& rhs);

 private:
  boost::dynamic_bitset<> flipped;
  boost::dynamic_bitset<> undecided;

  void flip(size_t rowcol);
  bool have_flipped(size_t row, size_t col) const;
};

#define RC2IND(row, col) col* nrows + row

template <class elem_t>
bool CartanEquiv::operator()(arma::Mat<elem_t> const& lhs,
                             arma::Mat<elem_t> const& rhs) {
  size_t const nrows = lhs.n_rows;
  // As we only check the lower triangle, assume that there are at most as many
  // columns as there are rows. In fact the matrices *should* be square.
  size_t const ncols = std::min(static_cast<size_t>(lhs.n_cols), nrows);
  bool result = nrows == rhs.n_rows && lhs.n_cols == rhs.n_cols;
  flipped.reset();
  flipped.resize(nrows, false);
  undecided.reset();
  undecided.resize(nrows, false);
  elem_t const* lptr = lhs.memptr();
  elem_t const* rptr = rhs.memptr();
  // Go down first column and determine which rows to flip. We decide that the
  // first column will be fixed, and so not flipped. We use this to determine
  // which columns do need flipping.
  // For matrices without any zeros this will 'fix' any sign discrepancies
  // between the matrices if they are equivalent. However any zeros in the first
  // column make this more difficult as 0 == -0 and 0 == 0.
  size_t col = 0;
  for (size_t row = 1; result && row < nrows; ++row) {
    elem_t const lval = lptr[RC2IND(row, col)];
    elem_t const rval = rptr[RC2IND(row, col)];
    if (lval == 0 && rval == 0) {
      undecided.set(row);
    } else if (lval == -rval) {
      // Flip the signs in the row
      flip(row);
    } else {
      result = (lval == rval);
    }
  }

  // Now need to deal with any zeros. Rows with a zero in the first column
  // cannot be determined to be flipped or not in the first pass above, so need
  // to travel down the column starting with a zero to see if all entries in
  // that column can be flipped.
  //
  // Check each undecided column to determine whether it should be flipped. All
  // previous columns have been fixed already, so we traverse the lower triangle
  // of the matrices.
  for (col = 1; result && col < ncols; ++col) {
    bool could_flip_col = undecided.test(col);
    bool need_flip_col = false;
    std::vector<size_t> simul_flip;
    std::vector<size_t> rows_to_flip;
    for (size_t row = col + 1; row < nrows; ++row) {
      elem_t const lval = lptr[RC2IND(row, col)];
      if (lval == 0) {
        continue;
      }
      elem_t const rval = rptr[RC2IND(row, col)];
      if (rval == 0) {
        result = false;
        break;
      }
      bool already_flipped = have_flipped(row, col);
      bool equal = (lval == rval);
      bool switched = (lval == -rval);
      result = (equal || switched);
      bool sign_matches =
          (already_flipped && switched) || (!already_flipped && equal);
      if (!sign_matches) {
        // Then need to flip either the row or the column
        if (undecided.test(row)) {
          rows_to_flip.push_back(row);
        } else {
          // Can't flip row, so need to flip column
          need_flip_col = true;
        }
      } else if (undecided.test(row)) {
        // Don't need to flip, but if we do flip the column, then we also need
        // to flip the row
        simul_flip.push_back(row);
      } else {
        // Shouldn't flip the column
        could_flip_col = false;
      }
    }
    if (could_flip_col) {
      flip(col);
      // Then also flip everything in simul_flip
      for (auto to_flip : simul_flip) {
        flip(to_flip);
        undecided.reset(to_flip);
      }
      // As we flipped col, cannot also flip the rows in rows_to_flip
      for (auto row : rows_to_flip) {
        undecided.reset(row);
      }
    } else {
      if (need_flip_col) {
        // Needed to flip col, but couldn't
        result = false;
      } else {
        // Didn't flip the column, so cannot flip the simul flip rows
        for (auto to_flip : simul_flip) {
          undecided.reset(to_flip);
        }
        // As we don't flip the column, we know which rows we *do* need to flip
        for (auto to_flip : rows_to_flip) {
          flip(to_flip);
          undecided.reset(to_flip);
        }
      }
    }
    // Note: don't need to call undecided.reset(col) as we will now move past it
    // and any later undecided rows just assume that previous ones have been
    // handled.
  }
  return result;
}
void CartanEquiv::flip(size_t rowcol) {
  flipped.flip(rowcol);
}
bool CartanEquiv::have_flipped(size_t row, size_t col) const {
  bool rval = flipped.test(row);
  bool cval = flipped.test(col);
  return rval != cval;
}
#undef RC2IND
}
#endif
