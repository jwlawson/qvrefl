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
	EXPECT_EQ(a, m.qv(0));
	EXPECT_EQ(a, m.qv(1));

	arma::mat b = { { 0, -1 }, { 1, 0 } };
	EXPECT_TRUE(util::equal(b, m.arma(0)));
	EXPECT_TRUE(util::equal(b, m.arma(1)));
}
TEST(MutStar, A3) {
	cluster::QuiverMatrix q("{ { 0 1 0 } { -1 0 1 } { 0 -1 0 } }");
	MutationStar m(q);

	cluster::QuiverMatrix a("{ { 0 -1 0 } { 1 0 1 } { 0 -1 0 } }");
	EXPECT_EQ(a, m.qv(0));
	cluster::QuiverMatrix b("{ { 0 -1 1 } { 1 0 -1 } { -1 1 0 } }");
	EXPECT_EQ(b, m.qv(1));
	cluster::QuiverMatrix c("{ { 0 1 0 } { -1 0 -1 } { 0 1 0 } }");
	EXPECT_EQ(c, m.qv(2));

	arma::mat d = { { 0, -1, 0 }, { 1, 0, 1 }, { 0, -1, 0 } };
	EXPECT_TRUE(util::equal(d, m.arma(0)));
	arma::mat e = { { 0, -1, 1 }, { 1, 0, -1 }, { -1, -1, 0 } };
	EXPECT_TRUE(util::equal(e, m.arma(1)));
	arma::mat f = { { 0, 1, 0 }, { -1, 0, -1 }, { 0, 1, 0 } };
	EXPECT_TRUE(util::equal(f, m.arma(2)));
}
}
