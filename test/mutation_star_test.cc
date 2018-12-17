/**
 * mutation_star_test.cc
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

#include "gtest/gtest.h"

#include "util.h"

namespace refl {
TEST(MutStar, A2) {
  cluster::QuiverMatrix q("{ { 0 1 } { -1 0 } }");
  MutationStar m(q);

  cluster::QuiverMatrix a("{ { 0 -1 } { 1 0 } }");
  EXPECT_TRUE(a.equals(m.qv(0)));
  EXPECT_TRUE(a.equals(m.qv(1)));

  arma::Mat<int> b = {{0, -1}, {1, 0}};
  EXPECT_TRUE(util::equal(b, m.arma(0)));
  EXPECT_TRUE(util::equal(b, m.arma(1)));
}
TEST(MutStar, A3) {
  cluster::QuiverMatrix q("{ { 0 1 0 } { -1 0 1 } { 0 -1 0 } }");
  MutationStar m(q);

  cluster::QuiverMatrix a("{ { 0 -1 0 } { 1 0 1 } { 0 -1 0 } }");
  EXPECT_TRUE(a.equals(m.qv(0)));
  cluster::QuiverMatrix b("{ { 0 -1 1 } { 1 0 -1 } { -1 1 0 } }");
  EXPECT_TRUE(b.equals(m.qv(1)));
  cluster::QuiverMatrix c("{ { 0 1 0 } { -1 0 -1 } { 0 1 0 } }");
  EXPECT_TRUE(c.equals(m.qv(2)));

  arma::Mat<int> d = {{0, -1, 0}, {1, 0, 1}, {0, -1, 0}};
  EXPECT_TRUE(util::equal(d, m.arma(0)));
  arma::Mat<int> e = {{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}};
  EXPECT_TRUE(util::equal(e, m.arma(1)));
  arma::Mat<int> f = {{0, 1, 0}, {-1, 0, -1}, {0, 1, 0}};
  EXPECT_TRUE(util::equal(f, m.arma(2)));
}
}  // namespace refl
