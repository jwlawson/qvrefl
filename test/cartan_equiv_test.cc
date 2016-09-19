/**
 * cartan_equiv_test.cc
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
#include "cartan_equiv.h"

#include "gtest/gtest.h"

namespace refl {
TEST(CartanEquiv, Same) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  CartanEquiv equiv;

  EXPECT_TRUE(equiv(a, a));

  arma::Mat<int> b{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(b, a));
}
TEST(CartanEquiv, OneSwitch) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{
      {2, -1, -1, -1}, {-1, 2, 1, 1}, {-1, 1, 2, 1}, {-1, 1, 1, 2}};
  arma::Mat<int> c{
      {2, -1, 1, 1}, {-1, 2, -1, -1}, {1, -1, 2, 1}, {1, -1, 1, 2}};
  arma::Mat<int> d{
      {2, 1, -1, 1}, {1, 2, 1, -1}, {-1, -1, 2, -1}, {1, 1, -1, 2}};
  arma::Mat<int> e{
      {2, 1, 1, -1}, {1, 2, 1, -1}, {1, 1, 2, -1}, {-1, -1, -1, 2}};

  CartanEquiv equiv;
  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(a, c));
  EXPECT_TRUE(equiv(a, d));
  EXPECT_TRUE(equiv(a, e));

  EXPECT_TRUE(equiv(b, a));
  EXPECT_TRUE(equiv(c, a));
  EXPECT_TRUE(equiv(d, a));
  EXPECT_TRUE(equiv(e, a));
}
TEST(CartanEquiv, TwoSwitch) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{
      {2, -1, 1, -1}, {-1, 2, -1, 1}, {1, -1, 2, -1}, {-1, 1, -1, 2}};
  arma::Mat<int> c{
      {2, -1, -1, 1}, {-1, 2, 1, -1}, {-1, 1, 2, -1}, {1, -1, -1, 2}};
  arma::Mat<int> d{
      {2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};

  CartanEquiv equiv;
  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(a, c));
  EXPECT_TRUE(equiv(a, d));

  EXPECT_TRUE(equiv(b, a));
  EXPECT_TRUE(equiv(c, a));
  EXPECT_TRUE(equiv(d, a));
}
TEST(CartanEquiv, NotSignEquiv) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{{2, 1, -1, -1}, {1, 2, 1, 1}, {-1, 1, 2, 1}, {-1, 1, 1, 2}};
  arma::Mat<int> c{
      {2, -1, 1, -1}, {-1, 2, -1, -1}, {1, -1, 2, 1}, {-1, -1, 1, 2}};
  arma::Mat<int> d{{2, 1, -1, 1}, {1, 2, 1, 1}, {-1, 1, 2, -1}, {1, 1, -1, 2}};
  arma::Mat<int> e{{2, 1, 1, 1}, {1, 2, 1, -1}, {1, 1, 2, -1}, {1, -1, -1, 2}};

  CartanEquiv equiv;
  EXPECT_FALSE(equiv(a, b));
  EXPECT_FALSE(equiv(a, c));
  EXPECT_FALSE(equiv(a, d));
  EXPECT_FALSE(equiv(a, e));

  EXPECT_FALSE(equiv(b, a));
  EXPECT_FALSE(equiv(c, a));
  EXPECT_FALSE(equiv(d, a));
  EXPECT_FALSE(equiv(e, a));
}
TEST(CartanEquiv, TrickySignNotEquiv) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{
      {2, 1, -1, -1}, {1, 2, 1, -1}, {-1, 1, 2, 1}, {-1, -1, 1, 2}};

  CartanEquiv equiv;
  EXPECT_FALSE(equiv(a, b));
  EXPECT_FALSE(equiv(b, a));
}
TEST(CartanEquiv, NotEqual) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{{2, 4, -3, 1}, {-3, 2, 2, 1}, {-1, 1, 2, 1}, {-1, 1, 1, 2}};
  CartanEquiv equiv;
  EXPECT_FALSE(equiv(a, b));
  EXPECT_FALSE(equiv(b, a));
}
TEST(CartanEquiv, DifferentSize) {
  arma::Mat<int> a{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
  arma::Mat<int> b{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  CartanEquiv equiv;
  EXPECT_FALSE(equiv(a, b));
  EXPECT_FALSE(equiv(b, a));
}
TEST(CartanEquiv, 5DimWithZeros) {
  arma::Mat<int> a{{2, 1, -1, -2, 0},
                   {1, 2, 1, -1, -1},
                   {-1, 1, 2, 1, -1},
                   {-2, -1, 1, 2, 0},
                   {0, -1, -1, 0, 2}};
  arma::Mat<int> b{{2, 1, 1, -2, 0},
                   {1, 2, -1, -1, -1},
                   {1, -1, 2, -1, 1},
                   {-2, -1, -1, 2, 0},
                   {0, -1, 1, 0, 2}};
  CartanEquiv equiv;

  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(b, a));
}
TEST(CartanEquiv, 5DimMoreZeros) {
  arma::Mat<int> a{{2, -2, 1, -1, 0},
                   {-2, 2, -1, 1, 0},
                   {1, -1, 2, 1, 1},
                   {-1, 1, 1, 2, 1},
                   {0, 0, 1, 1, 2}};
  arma::Mat<int> b{{2, 2, -1, 1, 0},
                   {2, 2, -1, 1, 0},
                   {-1, -1, 2, 1, 1},
                   {1, 1, 1, 2, 1},
                   {0, 0, 1, 1, 2}};
  CartanEquiv equiv;

  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(b, a));
}
TEST(CartanEquiv, MoreZeros) {
  arma::Mat<int> a{{2, 0, -1, 1, 0},
                   {0, 2, 1, -1, 0},
                   {-1, 1, 2, 0, 1},
                   {1, -1, 0, 2, 1},
                   {0, 0, 1, 1, 2}};

  // Multiply first row/col by -1
  arma::Mat<int> b{{2, 0, 1, -1, 0},
                   {0, 2, 1, -1, 0},
                   {1, 1, 2, 0, 1},
                   {-1, -1, 0, 2, 1},
                   {0, 0, 1, 1, 2}};
  CartanEquiv equiv;

  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(b, a));
}
TEST(CartanEquiv, YetMoreZeros) {
  arma::Mat<int> a{{2, 0, 1, 1, 0},
                   {0, 2, 1, 1, 0},
                   {1, 1, 2, 2, 1},
                   {1, 1, 2, 2, 1},
                   {0, 0, 1, 1, 2}};

  // Mutliply second row/col by -1
  arma::Mat<int> b{{2, 0, 1, 1, 0},
                   {0, 2, -1, -1, 0},
                   {1, -1, 2, 2, 1},
                   {1, -1, 2, 2, 1},
                   {0, 0, 1, 1, 2}};
  CartanEquiv equiv;

  EXPECT_TRUE(equiv(a, b));
  EXPECT_TRUE(equiv(b, a));
}
}
