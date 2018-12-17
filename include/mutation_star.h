/**
 * mutation_star.h
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
#ifndef REFL_MUTATION_STAR_H__
#define REFL_MUTATION_STAR_H__

#include "qv/quiver_matrix.h"

#include <armadillo>

namespace refl {
/**
 * Container of all quivers obtainable from the provided quiver by a single
 * mutation.
 *
 * These quivers can be accessed as either QuiverMatrix or arma::mat objects.
 */
class MutationStar {
 public:
  /**
   * Construct a new instance with the provided initial quiver.
   */
  MutationStar(cluster::QuiverMatrix const& q);
  /**
   * Get the mutation in the k-th direction of the initial quiver.
   */
  cluster::QuiverMatrix const& qv(uint_fast16_t k) const;
  arma::Mat<int> const& arma(uint_fast16_t k) const;

 private:
  std::vector<cluster::QuiverMatrix> _qv_vector;
  std::vector<arma::Mat<int>> _mat_vector;
};
}  // namespace refl
#endif
