/**
 * semi_positive_filter_test.cc
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
#include "semi_positive_filter.h"

#include "gtest/gtest.h"

namespace refl {
TEST(SemiPos, Id) {
  arma::Mat<int> a(5, 5, arma::fill::eye);
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, Not) {
  arma::Mat<int> a{{1, 2}, {2, 1}};
  SemiPositiveFilter f;
  EXPECT_FALSE(f(a));
}
TEST(SemiPos, ZeroEValue) {
  arma::Mat<int> a{{1, 1}, {1, 1}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, AllPos) {
  arma::Mat<int> a{{8, 5}, {5, 8}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, NegValue) {
  arma::Mat<int> a{{2, -1}, {-1, 2}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, Zeros) {
  arma::Mat<int> a{{0, 0}, {0, 0}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, DoubleEValuesWithNegative) {
  arma::Mat<int> a{{2, 5, -1}, {5, 2, 9}, {-1, 9, 2}};
  SemiPositiveFilter f;
  EXPECT_FALSE(f(a));
}
TEST(SemiPos, 3DimWithZero) {
  arma::Mat<int> a{{2, 1, -1}, {1, 2, 1}, {-1, 1, 2}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
TEST(SemiPos, DoubleAllPos) {
  arma::Mat<int> a{{2, 0, -1}, {0, 2, 1}, {-1, 1, 2}};
  SemiPositiveFilter f;
  EXPECT_TRUE(f(a));
}
}  // namespace refl
