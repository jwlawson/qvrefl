/**
 * cartan_iterator.h
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
#ifndef REFL_CARTAN_ITERATOR_H__
#define REFL_CARTAN_ITERATOR_H__

#include <armadillo>
#include <boost/dynamic_bitset.hpp>

#include "qv/quiver_matrix.h"

namespace refl {
/**
 * Iterator over all possible signed quasi-Cartan matrices for a given quiver.
 */
class CartanIterator {
 public:
  /**
   * Make  new instance which provides quasi-Cartan matrices for the given
   * quiver.
   */
  CartanIterator(cluster::QuiverMatrix const& q);
  /**
   * Check whether the iterator will return a valid matrix on the next call to
   * next()
   */
  bool has_next();
  /**
   * Get the next Cartan matrix.
   */
  arma::Mat<int>& next();

 private:
  uint_fast16_t const _number_vars;
  arma::Mat<int> const _initial;
  arma::Mat<int> _result;
  boost::dynamic_bitset<> const _zero_mask;
  uint_fast16_t const _non_zero_vars;
  uint_fast64_t _current_val;
  uint_fast64_t const _max_val;
};
}  // namespace refl
#endif
