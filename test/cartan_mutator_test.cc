/**
 * cartan_mutator_test.cc
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
#include "cartan_mutator.h"

#include "gtest/gtest.h"

namespace refl {
TEST(CartanMutator, A3) {
  cluster::QuiverMatrix a3{"{ { 0 1 -1 } { -1 0 1 } { 1 -1 0 } }"};
  arma::Mat<int> a{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
  CartanMutator mut(a3);

  arma::Mat<int> output;

  arma::Mat<int> exp1{{2, -1, 1}, {-1, 2, 0}, {1, 0, 2}};
  mut(a, 0, output);
  EXPECT_TRUE(std::equal(exp1.begin(), exp1.end(), output.begin()));

  arma::Mat<int> exp2{{2, 1, 0}, {1, 2, -1}, {0, -1, 2}};
  mut(a, 1, output);
  EXPECT_TRUE(std::equal(exp2.begin(), exp2.end(), output.begin()));

  arma::Mat<int> exp3{{2, 0, -1}, {0, 2, 1}, {-1, 1, 2}};
  mut(a, 2, output);
  EXPECT_TRUE(std::equal(exp3.begin(), exp3.end(), output.begin()));
}
}  // namespace refl
