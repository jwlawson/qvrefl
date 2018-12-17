#include "fully_compatible_check.h"
#include "gtest/gtest.h"

namespace refl {
TEST(FullyCompatible, TypeA) {
  cluster::QuiverMatrix m(
      "{ { 0 1 0 0 0 } "
      "  { -1 0 1 0 0 }"
      "	{ 0 -1 0 1 0 }"
      "	{ 0 0 -1 0 1 }"
      "	{ 0 0 0 -1 0 } }");
  arma::Mat<int> a{{2, 1, 0, 0, 0},
                   {1, 2, 1, 0, 0},
                   {0, 1, 2, 1, 0},
                   {0, 0, 1, 2, 1},
                   {0, 0, 0, 1, 2}};
  FullyCompatibleCheck chk;
  EXPECT_TRUE(chk(m, a));

  cluster::QuiverMatrix n(
      "{ { 0 1 0 0 0 } "
      "  { -1 0 1 0 0 }"
      " 	{ 0 -1 0 1 0 }"
      " 	{ 0 0 -1 0 1 }"
      " 	{ 0 0 0 -1 0 } }");
  arma::Mat<int> b{{2, -1, 0, 0, 0},
                   {-1, 2, -1, 0, 0},
                   {0, -1, 2, -1, 0},
                   {0, 0, -1, 2, -1},
                   {0, 0, 0, -1, 2}};
  EXPECT_TRUE(chk(n, b));
}
TEST(FullyCompatible, NotCompatible) {
  cluster::QuiverMatrix m(
      "{ { 0 1 -1 0 1 -1 }"
      " { -1 0 -1 1 0 1 }"
      " { 1 1 0 -1 -1 0 }"
      " { -1 0 1 -1 0 1 }"
      " { 1 -1 0 1 -1 0 } }");
  arma::Mat<int> a{{2, 1, 1, 0, 1, 1},
                   {1, 2, 1, 1, 0, 1},
                   {1, 1, 2, 1, 1, 0},
                   {1, 0, 1, 1, 2, 1},
                   {1, 1, 0, 1, 1, 2}};
  FullyCompatibleCheck chk;
  EXPECT_FALSE(chk(m, a));
}
TEST(FullyCompatible, CompatibleNotFully) {
  cluster::QuiverMatrix m(
      "{ { 0 1 1 0 -1 -1 } "
      "{ -1 0 -1 1 0 1 } "
      "{ -1 1 0 -1 1 0 } "
      "{ 0 -1 1 0 1 -1 } "
      "{ 1 0 -1 -1 0 1 } "
      "{ 1 -1 0 1 -1 0 } }");
  arma::Mat<int> a{{2, 1, 1, 0, -1, -1}, {-1, 2, -1, 1, 0, 1},
                   {-1, 1, 2, -1, 1, 0}, {0, -1, 1, 2, 1, -1},
                   {1, 0, -1, -1, 2, 1}, {1, -1, 0, 1, -1, 2}};
  FullyCompatibleCheck chk;
  EXPECT_FALSE(chk(m, a));
}
}  // namespace refl
