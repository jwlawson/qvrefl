/**
 * util_test.cc
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
#include "util.h"

#include "gtest/gtest.h"

namespace refl {
TEST(Equal, Small) {
  arma::Mat<int> a = {{1, 2, 3}, {2, 3, 1}, {9, 8, 7}};
  arma::Mat<int> b = {{1, 2, 3}, {2, 3, 1}, {9, 8, 7}};
  EXPECT_TRUE(util::equal(a, b));
}
TEST(Equal, NotSmall) {
  arma::Mat<int> a = {{1, 2, 3}, {4, 5, 6}, {9, 8, 7}};
  arma::Mat<int> b = {{1, 2, 3}, {4, 5, 5}, {9, 8, 7}};
  EXPECT_FALSE(util::equal(a, b));
}
TEST(Equal, Close) {
  arma::Mat<int> a = {{0, 2}, {4, 5}};
  arma::Mat<int> b = {{0, 1}, {4, 5}};
  EXPECT_FALSE(util::equal(a, b));
}
TEST(Equal, Sizes) {
  arma::Mat<int> a = {{0, 1}, {4, 5}};
  arma::Mat<int> b = {{0, 1, 0}, {4, 5, 0}, {0, 0, 0}};
  EXPECT_FALSE(util::equal(a, b));
}
TEST(Adaptor, 2Dim) {
  cluster::QuiverMatrix q("{ { 0 2 } { -2 0 } }");
  arma::Mat<int> exp = {{0, 2}, {-2, 0}};
  arma::Mat<int> res = util::to_arma(q);
  EXPECT_TRUE(util::equal(exp, res));
}
TEST(Adaptor, 3Dim) {
  cluster::QuiverMatrix q("{ { 0 2 0 } { -2 0 1 } { 0 -1 0 } }");
  arma::Mat<int> exp = {{0, 2, 0}, {-2, 0, 1}, {0, -1, 0}};
  arma::Mat<int> res = util::to_arma(q);
  EXPECT_TRUE(util::equal(exp, res));
}
TEST(Gram, IdId) {
  arma::Mat<int> v = {{1, 0}, {0, 1}};
  arma::Mat<int> res = util::gram(v, v);
  EXPECT_TRUE(util::equal(v, res));
}
TEST(Gram, Id) {
  arma::Mat<int> v = {{2, 1}, {1, 1}};
  arma::Mat<int> prod = {{1, 0}, {0, 1}};
  arma::Mat<int> exp = {{5, 3}, {3, 2}};
  arma::Mat<int> res = util::gram(v, prod);
  EXPECT_TRUE(util::equal(exp, res));
}
TEST(Gram, Prod) {
  arma::Mat<int> v = {{2, 1}, {1, 1}};
  arma::Mat<int> prod = {{-1, 2}, {2, 1}};
  arma::Mat<int> exp = {{5, 5}, {5, 4}};
  arma::Mat<int> res = util::gram(v, prod);
  EXPECT_TRUE(util::equal(exp, res));
}
}  // namespace refl
