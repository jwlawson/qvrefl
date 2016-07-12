/**
 * unique_matrix_filter_test.cc
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
	arma::Mat<int> b(4, 4,  arma::fill::eye);
	EXPECT_TRUE(f(b));
}
TEST(UniqueMat, SameVals) {
	arma::Mat<int> a { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	UniqueMatrixFilter f;
	EXPECT_TRUE(f(a));
	arma::Mat<int> b { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	EXPECT_FALSE(f(b));
}
TEST(UniqueMat, DiffSigns) {
	arma::Mat<int> a { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	UniqueMatrixFilter f;
	EXPECT_TRUE(f(a));
	arma::Mat<int> b { { 1, 2, 3 }, { 4, 5, -6 }, { 7, 8, 9 } };
	EXPECT_TRUE(f(b));
}
TEST(UniqueMat, DiffVals) {
	arma::Mat<int> a { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	UniqueMatrixFilter f;
	EXPECT_TRUE(f(a));
	arma::Mat<int> b { { 1, 2, 1 }, { 1, 5, 6 }, { 7, 8, 9 } };
	EXPECT_TRUE(f(b));
}
TEST(UniqueMat, DiffLastVals) {
	arma::Mat<int> a { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	UniqueMatrixFilter f;
	EXPECT_TRUE(f(a));
	arma::Mat<int> b { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 1 } };
	EXPECT_TRUE(f(b));
}
TEST(UniqueMat, Permutation) {
	arma::Mat<int> a { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	UniqueMatrixFilter f;
	EXPECT_TRUE(f(a));
	arma::Mat<int> b { { 5, 4, 6 }, { 2, 1, 3 }, { 8, 7, 9 } };
	EXPECT_TRUE(f(b));
}
}
