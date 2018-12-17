/**
 * unique_matrix_filter_test.cc
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
#include "gtest/gtest.h"

namespace refl {
TEST(UniqueMat, SameMat) {
  arma::Mat<int> a(3, 3, arma::fill::eye);
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  EXPECT_FALSE(f(a));
  EXPECT_FALSE(f(a));
  EXPECT_FALSE(f(a));
}
TEST(UniqueMat, DiffSizes) {
  arma::Mat<int> a(3, 3, arma::fill::eye);
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b(4, 4, arma::fill::eye);
  EXPECT_TRUE(f(b));
}
TEST(UniqueMat, SameVals) {
  arma::Mat<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  EXPECT_FALSE(f(b));
}
TEST(UniqueMat, DiffSigns) {
  arma::Mat<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b{{1, 2, 3}, {4, 5, -6}, {7, 8, 9}};
  EXPECT_TRUE(f(b));
}
TEST(UniqueMat, DiffVals) {
  arma::Mat<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b{{1, 2, 1}, {1, 5, 6}, {7, 8, 9}};
  EXPECT_TRUE(f(b));
}
TEST(UniqueMat, DiffLastVals) {
  arma::Mat<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b{{1, 2, 3}, {4, 5, 6}, {7, 8, 1}};
  EXPECT_TRUE(f(b));
}
TEST(UniqueMat, Permutation) {
  arma::Mat<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  UniqueMatrixFilter f;
  EXPECT_TRUE(f(a));
  arma::Mat<int> b{{5, 4, 6}, {2, 1, 3}, {8, 7, 9}};
  EXPECT_TRUE(f(b));
}
}  // namespace refl
