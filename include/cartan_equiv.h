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
class CartanEquiv {
public:
  template <class elem_t>
  bool operator()(arma::Mat<elem_t> const &lhs, arma::Mat<elem_t> const &rhs);

private:
  boost::dynamic_bitset<> flipped;
  boost::dynamic_bitset<> undecided;

  void flip(size_t rowcol, size_t nrows, size_t ncols);
};

#define RC2IND(row, col) col *nrows + row

template <class elem_t>
bool CartanEquiv::operator()(arma::Mat<elem_t> const &lhs,
                             arma::Mat<elem_t> const &rhs) {
  size_t const nrows = lhs.n_rows;
  size_t const ncols = lhs.n_cols;
  bool result = nrows == rhs.n_rows && ncols == rhs.n_cols;
  flipped.reset();
  flipped.resize(nrows * ncols, false);
  undecided.reset();
  undecided.resize(nrows, false);
  elem_t const *lptr = lhs.memptr();
  elem_t const *rptr = rhs.memptr();
  // Go down first column and determine which rows to flip.
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
      flip(row, nrows, ncols);
    } else {
      result = (lval == rval);
    }
  }

  // Now need to deal with any zeros. Rows with a zero in the first column
  // cannot be determined to be flipped or not in the first pass above, so need
  // to travel down the column starting with a zero to see if all entries in
  // that column can be flipped.
  //
  // We do this in two sections, the upper triangle and then the lower triangle.
  // For the upper triangle section we just check whether the values in the
  // matrix need flipping, whereas in the lower triangle we also need to be
  // aware that we could simultaneously flip both this current column and a
  // column further in the matrix, so these possible future flips are tracked in
  // simul_flip.
  for (col = 1; result && col < ncols; ++col) {
    bool could_flip = undecided.test(col);
    std::vector<size_t> simul_flip;
    for (size_t row = 1; could_flip && row < col; ++row) {
      elem_t const lval = lptr[RC2IND(row, col)];
      if (lval == 0) {
        continue;
      }
      elem_t const rval = rptr[RC2IND(row, col)];
      bool already_flipped = flipped.test(RC2IND(col, row));
      bool sign_matches = (already_flipped && lval == -rval) ||
                          (!already_flipped && lval == rval);
      if (sign_matches) {
        // When checking the upper triangle, there can be no previous rows to
        // flip (i.e. undecided.test(row) is always false).
        could_flip = false;
      }
    }
    for (size_t row = col + 1; could_flip && row < nrows; ++row) {
      elem_t const lval = lptr[RC2IND(row, col)];
      if (lval == 0) {
        continue;
      }
      elem_t const rval = rptr[RC2IND(row, col)];
      bool already_flipped = flipped.test(RC2IND(row, col));
      bool sign_matches = (already_flipped && lval == -rval) ||
                          (!already_flipped && lval == rval);
      if (!sign_matches) {
        // Then need to flip either the row or the column
      } else if (undecided.test(row)) {
        // Don't need to flip, but could if we flip both the row and column
        simul_flip.push_back(row);
      } else {
        // Don't need to flip
        could_flip = false;
      }
    }
    if (could_flip) {
      flip(col, nrows, ncols);
      // Then also flip everything in simul_flip
      for (auto sflip = simul_flip.begin(), end = simul_flip.end();
           sflip != end; ++sflip) {
        size_t to_flip = *sflip;
        flip(to_flip, nrows, ncols);
        undecided.reset(to_flip);
      }
    }
    // Note: don't need to call undecided.reset(col) as we will now move past it
    // and any later undecided rows just assume that previous ones have been
    // handled.
  }
  // At this point all flips have been decided, so run through the matrix and
  // check that everything matches up.
  //
  // We could include some of this checking in the above, so that the code can
  // fail faster, but that probably will cause it to be even more inpenetrable.
  for (col = 1; result && col < ncols; ++col) {
    for (size_t row = col + 1; result && row < nrows; ++row) {
      elem_t const lval = lptr[RC2IND(row, col)];
      elem_t const rval = rptr[RC2IND(row, col)];
      if (lval == 0 && rval == 0) {
        continue;
      }
      bool already_flipped = flipped.test(RC2IND(row, col));
      result = (already_flipped && lval == -rval) ||
               (!already_flipped && lval == rval);
    }
  }

  return result;
}
void CartanEquiv::flip(size_t rowcol, size_t nrows, size_t ncols) {
  for (size_t second = 1; second < rowcol; ++second) {
    flipped.flip(RC2IND(rowcol, second));
  }
  for (size_t first = rowcol + 1; first < ncols; ++first) {
    flipped.flip(RC2IND(first, rowcol));
  }
}
#undef RC2IND
}
#endif
