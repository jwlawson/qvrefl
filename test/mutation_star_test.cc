/**
 * mutation_star_test.cc
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

	arma::Mat<int> b = { { 0, -1 }, { 1, 0 } };
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

	arma::Mat<int> d = { { 0, -1, 0 }, { 1, 0, 1 }, { 0, -1, 0 } };
	EXPECT_TRUE(util::equal(d, m.arma(0)));
	arma::Mat<int> e = { { 0, -1, 1 }, { 1, 0, -1 }, { -1, -1, 0 } };
	EXPECT_TRUE(util::equal(e, m.arma(1)));
	arma::Mat<int> f = { { 0, 1, 0 }, { -1, 0, -1 }, { 0, 1, 0 } };
	EXPECT_TRUE(util::equal(f, m.arma(2)));
}
}
