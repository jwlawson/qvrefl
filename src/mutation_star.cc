/**
 * mutation_star.cc
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
#include "mutation_star.h"

#include "util.h"

namespace refl {
MutationStar::MutationStar(cluster::QuiverMatrix const& q)
    : _qv_vector(q.num_cols(),
                 cluster::QuiverMatrix(q.num_rows(), q.num_cols()))
    , _mat_vector() {
  _mat_vector.reserve(q.num_cols());
  for (int_fast16_t i = 0; i < q.num_cols(); ++i) {
    q.mutate(i, _qv_vector[i]);
    _mat_vector.push_back(util::to_arma(_qv_vector[i]));
  }
}

cluster::QuiverMatrix const&
MutationStar::qv(uint_fast16_t k) const {
  return _qv_vector[k];
}
arma::Mat<int> const&
MutationStar::arma(uint_fast16_t k) const {
  return _mat_vector[k];
}
}  // namespace refl
