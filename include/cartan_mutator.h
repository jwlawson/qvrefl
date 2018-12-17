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
  CartanMutator(cluster::QuiverMatrix const &quiver);

  template <class elem_t>
  void operator()(arma::Mat<elem_t> const &cartan, uint_fast16_t k,
                  arma::Mat<elem_t> &output);

 private:
  cluster::QuiverMatrix const &m_quiver;
};

inline CartanMutator::CartanMutator(cluster::QuiverMatrix const &q)
    : m_quiver(q) {}

template <class elem_t>
void CartanMutator::operator()(arma::Mat<elem_t> const &cartan, uint_fast16_t k,
                               arma::Mat<elem_t> &output) {
#ifdef __cpp_constexpr
#if __cpp_constexpr >= 201603
#define CONSTEXPR_LAMBDA constexpr
#else
#define CONSTEXPR_LAMBDA
#endif
#endif
  static CONSTEXPR_LAMBDA auto sgn = [](auto const &x) {
    return x < 0 ? -1 : 1;
  };
  static CONSTEXPR_LAMBDA auto min0 = [](auto const &x) {
    return x < 0 ? 0 : x;
  };
  output.set_size(arma::size(cartan));
  uint_fast16_t max = cartan.n_rows;
  for (uint_fast16_t col = 0; col < max; ++col) {
    uint_fast16_t row = 0;
    for (; row < col; ++row) {
      output.at(row, col) = output.at(col, row);
    }
    output.at(row, col) = 2;
    ++row;
    for (; row < max; ++row) {
      if (col == k) {
        output.at(row, col) = sgn(m_quiver.get(row, k)) * cartan.at(row, k);
      } else if (row == k) {
        output.at(row, col) = -sgn(m_quiver.get(k, col)) * cartan.at(k, col);
      } else {
        output.at(row, col) =
            cartan.at(row, col) -
            sgn(cartan.at(row, k) * cartan.at(k, col)) *
                min0(m_quiver.get(row, k) * m_quiver.get(k, col));
      }
    }
  }
#undef CONSTEXPR_LAMBDA
}

}  // namespace refl

#endif  // REFL_CARTAN_MUTATOR_H__
