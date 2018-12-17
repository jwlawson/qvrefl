/**
 * cartan_iterator_test.cc
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
#include "cartan_iterator.h"

#include "gtest/gtest.h"

#include "util.h"

namespace refl {
TEST(CartanIt, 2Dim) {
  cluster::QuiverMatrix q("{ { 0 1 } { -1 0 } }");
  CartanIterator c(q);
  ASSERT_TRUE(c.has_next());
  arma::Mat<int>& a = c.next();
  arma::Mat<int> exp = {{2, 1}, {1, 2}};
  EXPECT_TRUE(util::equal(a, exp));
  ASSERT_TRUE(c.has_next());
  a = c.next();
  arma::Mat<int> exp2 = {{2, -1}, {-1, 2}};
  EXPECT_TRUE(util::equal(a, exp2));
  ASSERT_FALSE(c.has_next());
}
TEST(CartanIt, Count4) {
  cluster::QuiverMatrix q(
      "{ { 0 1 1 1 } { -1 0 1 1 } { -1 -1 0 1 } { -1 -1 -1 0 } }");
  CartanIterator c(q);
  ASSERT_TRUE(c.has_next());
  int count = 0;
  while (c.has_next()) {
    c.next();
    ++count;
  }
  EXPECT_EQ(64, count);
}
}  // namespace refl
