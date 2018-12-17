/**
 * fully_compatible_check.h
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
#ifndef REFL_FULLY_COMPATIBLE_CHECK_H__
#define REFL_FULLY_COMPATIBLE_CHECK_H__

#include "mutation_star.h"

#include "qv/quiver_matrix.h"

#include <armadillo>

namespace refl {
struct FullyCompatibleCheck {
  bool operator()(cluster::QuiverMatrix const& quiver,
                  arma::Mat<int> const& cartan);
  bool operator()(cluster::QuiverMatrix const& quiver, MutationStar const& star,
                  arma::Mat<int> const& cartan);
};
struct FixedQuiverFullyCompatibleCheck {
  FixedQuiverFullyCompatibleCheck(cluster::QuiverMatrix q);
  bool operator()(arma::Mat<int> const& cartan);

 private:
  cluster::QuiverMatrix quiver;
  MutationStar star;
  FullyCompatibleCheck check;
};
}  // namespace refl
#endif
