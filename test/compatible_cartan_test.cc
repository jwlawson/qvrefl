/**
 * compatible_cartan_test.cc
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
#include "compatible_cartan.h"

#include "gtest/gtest.h"

namespace refl {
TEST(CC, PosSame) {
  arma::mat a = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
  cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
  CompatibleCartan cc;
  EXPECT_TRUE(cc(q, a));
}
TEST(CC, NegSame) {
  arma::mat a = {{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}};
  cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
  CompatibleCartan cc;
  EXPECT_TRUE(cc(q, a));
}
TEST(CC, MixSame) {
  arma::mat a = {{2, -1, 1}, {-1, 2, 1}, {1, 1, 2}};
  cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
  CompatibleCartan cc;
  EXPECT_TRUE(cc(q, a));
}
TEST(CC, Diff1) {
  arma::mat a = {{2, -1, 1}, {-1, 2, 1}, {1, 1, 2}};
  cluster::QuiverMatrix q("{ { 0 2 1 } { -2 0 1 } { -1 -1 0 } }");
  CompatibleCartan cc;
  EXPECT_FALSE(cc(q, a));
}
TEST(CC, Diff2) {
  arma::mat a = {{2, -1, 1}, {-1, 2, 1}, {1, 1, 2}};
  cluster::QuiverMatrix q("{ { 0 1 2 } { -1 0 1 } { -2 -1 0 } }");
  CompatibleCartan cc;
  EXPECT_FALSE(cc(q, a));
}
}  // namespace refl
