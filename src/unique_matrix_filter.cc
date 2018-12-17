/**
 * unique_matrix_filter.cc
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
#include "unique_matrix_filter.h"

namespace refl {
bool
UniqueMatrixFilter::operator()(arma::Mat<int> const& m) {
  return _matrices.insert(m).second;
}
std::size_t
UniqueMatrixFilter::MatrixHash::operator()(arma::Mat<int> const& m) const {
  std::size_t result = 17;
  uint_fast64_t exp = 31;
  for (uint_fast64_t i = 0; i < m.n_elem; ++i) {
    result *= exp;
    result += m(i);
  }
  return result;
}
bool
UniqueMatrixFilter::MatrixEqual::operator()(arma::Mat<int> const& lhs,
                                            arma::Mat<int> const& rhs) const {
  bool result = lhs.n_cols == rhs.n_cols && lhs.n_rows == rhs.n_rows;
  for (uint_fast64_t i = 0; result && i < lhs.n_elem; ++i) {
    result = lhs(i) == rhs(i);
  }
  return result;
}
}  // namespace refl
