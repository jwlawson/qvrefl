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

namespace refl {
class CartanEquiv {
public:
  template <class elem_t>
  bool operator()(arma::Mat<elem_t> const &lhs, arma::Mat<elem_t> const &rhs);

private:
  boost::dynamic_bitset<> flipped;
};

template <class elem_t>
bool CartanEquiv::operator()(arma::Mat<elem_t> const &lhs,
                             arma::Mat<elem_t> const &rhs) {
  size_t const nrows = lhs.n_rows;
  size_t const ncols = lhs.n_cols;
  bool result = nrows == rhs.n_rows && ncols == rhs.n_cols;
  auto fl_ind = [nrows](size_t row, size_t col) { return col * nrows + row; };
  flipped.reset();
  flipped.resize(nrows * ncols, false);
  size_t col = 0;
  for (size_t row = col + 1; result && row < nrows; ++row) {
    bool already_flipped = flipped.test(fl_ind(row, col));
    if (lhs(row, col) == -rhs(row, col)) {
      // If already flipped, then this has been fixed.
      if (!already_flipped) {
        for (size_t second = col + 1; second < row; ++second) {
          flipped.flip(fl_ind(row, second));
        }
        for (size_t first = row + 1; first < ncols; ++first) {
          flipped.flip(fl_ind(first, row));
        }
      }
    } else {
      result = !already_flipped && lhs(row, col) == rhs(row, col);
    }
  }
  for (col = 1; result && col < ncols; ++col) {
    for (size_t row = col + 1; result && row < nrows; ++row) {
      bool already_flipped = flipped.test(fl_ind(row, col));
      result = (already_flipped && lhs(row, col) == -rhs(row, col)) ||
               (!already_flipped && lhs(row, col) == rhs(row, col));
    }
  }

  return result;
}
}
#endif
