/**
 * util_test.cc
 */
#include "util.h"

#include "gtest/gtest.h"

namespace refl {
TEST(Equal, Small) {
	arma::Mat<int> a = { { 1, 2, 3 }, { 2, 3, 1 }, { 9, 8, 7 } };
	arma::Mat<int> b = { { 1, 2, 3 }, { 2, 3, 1 }, { 9, 8, 7 } };
	EXPECT_TRUE(util::equal(a, b));
}
TEST(Equal, NotSmall) {
	arma::Mat<int> a = { { 1, 2, 3 }, { 4, 5, 6 }, { 9, 8, 7 } };
	arma::Mat<int> b = { { 1, 2, 3 }, { 4, 5, 5 }, { 9, 8, 7 } };
	EXPECT_FALSE(util::equal(a, b));
}
TEST(Equal, Close) {
	arma::Mat<int> a = { { 0, 2 }, { 4, 5 } };
	arma::Mat<int> b = { { 0, 1 }, { 4, 5 } };
	EXPECT_FALSE(util::equal(a, b));
}
TEST(Equal, Sizes) {
	arma::Mat<int> a = { { 0, 1 }, { 4, 5 } };
	arma::Mat<int> b = { { 0, 1, 0 }, { 4, 5, 0 }, { 0, 0, 0 } };
	EXPECT_FALSE(util::equal(a, b));
}
TEST(Adaptor, 2Dim) {
	cluster::QuiverMatrix q("{ { 0 2 } { -2 0 } }");
	arma::Mat<int> exp = { { 0, 2 }, { -2, 0 } };
	arma::Mat<int> res = util::to_arma(q);
	EXPECT_TRUE(util::equal(exp, res));
}
TEST(Adaptor, 3Dim) {
	cluster::QuiverMatrix q("{ { 0 2 0 } { -2 0 1 } { 0 -1 0 } }");
	arma::Mat<int> exp = { { 0, 2, 0 }, { -2, 0, 1 }, { 0, -1, 0 } };
	arma::Mat<int> res = util::to_arma(q);
	EXPECT_TRUE(util::equal(exp, res));
}
TEST(Gram, IdId) {
	arma::Mat<int> v = { { 1, 0 }, { 0, 1 } };
	arma::Mat<int> res = util::gram(v, v);
	EXPECT_TRUE(util::equal(v, res));
}
TEST(Gram, Id) {
	arma::Mat<int> v = { { 2, 1 }, { 1, 1 } };
	arma::Mat<int> prod = { { 1, 0 }, { 0, 1 } };
	arma::Mat<int> exp = { { 5, 3 }, { 3, 2 } };
	arma::Mat<int> res = util::gram(v, prod);
	EXPECT_TRUE(util::equal(exp, res));
}
TEST(Gram, Prod) {
	arma::Mat<int> v = { { 2, 1 }, { 1, 1 } };
	arma::Mat<int> prod = { { -1, 2 }, { 2, 1 } };
	arma::Mat<int> exp = { { 5, 5 }, { 5, 4 } };
	arma::Mat<int> res = util::gram(v, prod);
	EXPECT_TRUE(util::equal(exp, res));
}
}
