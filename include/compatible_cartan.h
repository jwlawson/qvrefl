/**
 * compatible_cartan.h
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
#ifndef REFL_COMPATIBLE_CARTAN_H__
#define REFL_COMPATIBLE_CARTAN_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
class CompatibleCartan {
 public:
  template <class mat>
  bool operator()(cluster::QuiverMatrix const& lhs, mat const& rhs) const;
};
template <class mat>
bool
CompatibleCartan::operator()(cluster::QuiverMatrix const& lhs,
                             mat const& rhs) const {
  bool result = static_cast<uint_fast16_t>(lhs.num_cols()) == rhs.n_cols &&
                static_cast<uint_fast16_t>(lhs.num_rows()) == rhs.n_rows;
  /* As the two types of matrix are arranged differently in memory, we cannot
   * use quick comparisons, but have to go through each value separately. */
  for (uint_fast16_t row = 0; result && row < rhs.n_rows; ++row) {
    for (uint_fast16_t col = 0; result && col < rhs.n_cols; ++col) {
      if (row == col) continue;
      result = std::abs(lhs.get(row, col)) == std::abs(rhs(row, col));
    }
  }
  return result;
}
}  // namespace refl
#endif
