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
  std::vector<size_t> undecided;
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
  undecided.clear();
  size_t col = 0;
  for (size_t row = col + 1; result && row < nrows; ++row) {
    elem_t const lval = lhs(row, col);
    elem_t const rval = rhs(row, col);
    if (lval == 0 && rval == 0) {
      undecided.push_back(row);
    } else if (lval == -rval) {
      for (size_t second = 1; second < row; ++second) {
        flipped.flip(fl_ind(row, second));
      }
      for (size_t first = row + 1; first < ncols; ++first) {
        flipped.flip(fl_ind(first, row));
      }
    } else {
      result = lval == rval;
    }
  }
  for (col = 1; result && col < ncols; ++col) {
    for (size_t row = col + 1; result && row < nrows; ++row) {
      elem_t const lval = lhs(row, col);
      elem_t const rval = rhs(row, col);
      if (lval == 0 && rval == 0) {
        continue;
      }
      bool already_flipped = flipped.test(fl_ind(row, col));
      result = (already_flipped && lval == -rval) ||
               (!already_flipped && lval == rval);
      auto iter_pos = std::find(undecided.begin(), undecided.end(), row);
			size_t flip_row = row;
      if (iter_pos == undecided.end()) {
        iter_pos = std::find(undecided.begin(), undecided.end(), col);
				flip_row = col;
      }
      if (iter_pos != undecided.end()) {
        if (!result) {
          for (size_t second = 1; second < flip_row; ++second) {
            flipped.flip(fl_ind(flip_row, second));
          }
          for (size_t first = flip_row + 1; first < ncols; ++first) {
            flipped.flip(fl_ind(first, flip_row));
          }
        }
        undecided.erase(iter_pos);
        result = true;
      }
    }
  }

  return result;
}
}
#endif
