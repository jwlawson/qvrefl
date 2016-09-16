/**
 * cartan_mutator.h
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
#ifndef REFL_CARTAN_MUTATOR_H__
#define REFL_CARTAN_MUTATOR_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
class CartanMutator {
 public:
  CartanMutator(cluster::QuiverMatrix const& quiver);

  template <class elem_t>
  void operator()(arma::Mat<elem_t> const& cartan,
                                 uint_fast16_t k, arma::Mat<elem_t>& output);

 private:
  cluster::QuiverMatrix const& _quiver;
};

CartanMutator::CartanMutator(cluster::QuiverMatrix const& q)
    : _quiver(q) {
}

template <class elem_t>
void CartanMutator::operator()(arma::Mat<elem_t> const& cartan, uint_fast16_t k,
                               arma::Mat<elem_t>& output) {
  static constexpr auto sgn = [](auto const& x) { return x < 0 ? -1 : 1; };
  static constexpr auto min0 = [](auto const& x) { return x < 0 ? 0 : x; };
	output.set_size(arma::size(cartan));
  for (uint_fast16_t row = 0, max = cartan.n_rows; row < max; ++row) {
    for (uint_fast16_t col = 0; col < max; ++col) {
      if (row == col) {
        output.at(row, col) = 2;
      } else if (col == k) {
        output.at(row, col) = sgn(_quiver.get(row, k)) * cartan.at(row, k);
      } else if (row == k) {
        output.at(row, col) = -sgn(_quiver.get(k, col)) * cartan.at(k, col);
      } else {
        output.at(row, col) =
            cartan.at(row, col) -
            sgn(cartan.at(row, k) * cartan.at(k, col)) *
                min0(_quiver.get(row, k) * _quiver.get(k, col));
      }
    }
  }
}
}
#endif
