/**
 * util_test.cc
 */
#include "util.h"

#include "gtest/gtest.h"

namespace refl {
namespace {
bool
equal(arma::mat const& a, arma::mat const& b) {
	constexpr double tol = 1e-10;
	bool result = a.n_cols == b.n_cols && a.n_rows == b.n_rows;
	for(auto a_it = a.begin(), a_end = a.end(), b_it = b.begin(), b_end = b.end();
			a_it != a_end && b_it != b_end;
			++a_it, ++b_it) {
		result = result && std::abs(*a_it - *b_it) < tol;
	}
	return result;
}
}

TEST(Adaptor, 2Dim) {
	cluster::QuiverMatrix q("{ { 0 2 } { -2 0 } }");
	arma::mat exp = { { 0, 2 }, { -2, 0 } };
	arma::mat res = util::to_arma(q);
	EXPECT_TRUE(equal(exp, res));
}
TEST(Adaptor, 3Dim) {
	cluster::QuiverMatrix q("{ { 0 2 0 } { -2 0 1 } { 0 -1 0 } }");
	arma::mat exp = { { 0, 2, 0 }, { -2, 0, 1 }, { 0, -1, 0 } };
	arma::mat res = util::to_arma(q);
	EXPECT_TRUE(equal(exp, res));
}
TEST(Gram, IdId) {
	arma::mat v = { { 1, 0 }, { 0, 1 } };
	arma::mat res = util::gram(v, v);
	EXPECT_TRUE(equal(v, res));
}
TEST(Gram, Id) {
	arma::mat v = { { 2, 1 }, { 1, 1 } };
	arma::mat prod = { { 1, 0 }, { 0, 1 } };
	arma::mat exp = { { 5, 3 }, { 3, 2 } };
	arma::mat res = util::gram(v, v);
	EXPECT_TRUE(equal(exp, res));
}
TEST(Gram, Prod) {
	arma::mat v = { { 2, 1 }, { 1, 1 } };
	arma::mat prod = { { -1, 2 }, { 2, 1 } };
	arma::mat exp = { { 5, 5 }, { 5, 4 } };
	arma::mat res = util::gram(v, v);
	EXPECT_TRUE(equal(exp, res));
}
}
