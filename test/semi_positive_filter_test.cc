/**
 * semi_positive_filter_test.cc
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
	arma::Mat<int> a { { 1, 2 }, { 2, 1 } };
	SemiPositiveFilter f;
	EXPECT_FALSE(f(a));
}
TEST(SemiPos, ZeroEValue) {
	arma::Mat<int> a { { 1, 1 }, { 1, 1 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
TEST(SemiPos, AllPos) {
	arma::Mat<int> a { { 8, 5 }, { 5, 8 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
TEST(SemiPos, NegValue) {
	arma::Mat<int> a { { 2, -1 }, { -1, 2 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
TEST(SemiPos, Zeros) {
	arma::Mat<int> a { { 0, 0 }, { 0, 0 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
TEST(SemiPos, DoubleEValuesWithNegative) {
	arma::Mat<int> a { { 2, 5, -1 }, { 5, 2, 9 }, { -1, 9, 2 } };
	SemiPositiveFilter f;
	EXPECT_FALSE(f(a));
}
TEST(SemiPos, 3DimWithZero) {
	arma::Mat<int> a { { 2, 1, -1 }, { 1, 2, 1 }, { -1, 1, 2 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
TEST(SemiPos, DoubleAllPos) {
	arma::Mat<int> a { { 2, 0, -1 }, { 0, 2, 1 }, { -1, 1, 2 } };
	SemiPositiveFilter f;
	EXPECT_TRUE(f(a));
}
}

