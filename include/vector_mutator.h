/**
 * vector_mutator.h
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
#ifndef REFL_VECTOR_MUTATOR_H__
#define REFL_VECTOR_MUTATOR_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
class VectorMutator {
public:
  /**
   * Construct a mutator with quasi-Cartan matrix A for quiver Q.
   */
  VectorMutator(cluster::QuiverMatrix const &Q, arma::Mat<int> const &A);
  /**
   * Mutate the vectors in the direction k, using the set quasi-Cartan matrix
   * and quiver. The resulting vectors are computed in the provided output
   * matrix.
   */
  template <class elem_t>
  bool mutate(arma::Mat<elem_t> const &vectors, uint_fast16_t k,
              arma::Mat<elem_t> &output);

private:
  cluster::QuiverMatrix const &m_quiver;
  arma::Mat<int> const &m_cartan;
};

VectorMutator::VectorMutator(cluster::QuiverMatrix const &Q,
                             arma::Mat<int> const &A)
    : m_quiver(Q), m_cartan(A) {}

template <class elem_t>
bool VectorMutator::mutate(arma::Mat<elem_t> const &vectors, uint_fast16_t k,
                           arma::Mat<elem_t> &output) {
  for (int_fast16_t i = 0; i < m_quiver.num_rows(); ++i) {
    if (m_quiver.get(i, k) > 0) {
      output.col(i) = vectors.col(i) - m_cartan(i, k) * vectors.col(k);
    } else {
      output.col(i) = vectors.col(i);
    }
  }
  return true;
}
}
#endif
