/**
 * compatible_cartan_test.cc
 */
#include "compatible_cartan.h"

#include "gtest/gtest.h"

namespace refl {
TEST(CC, PosSame) {
	arma::mat a = { { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
	CompatibleCartan cc;
	EXPECT_TRUE(cc(q, a));
}
TEST(CC, NegSame) {
	arma::mat a = { { 2, -1, -1 }, { -1, 2, -1 }, { -1, -1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
	CompatibleCartan cc;
	EXPECT_TRUE(cc(q, a));
}
TEST(CC, MixSame) {
	arma::mat a = { { 2, -1, 1 }, { -1, 2, 1 }, { 1, 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 1 } { -1 0 1 } { -1 -1 0 } }");
	CompatibleCartan cc;
	EXPECT_TRUE(cc(q, a));
}
TEST(CC, Diff1) {
	arma::mat a = { { 2, -1, 1 }, { -1, 2, 1 }, { 1, 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 2 1 } { -2 0 1 } { -1 -1 0 } }");
	CompatibleCartan cc;
	EXPECT_FALSE(cc(q, a));
}
TEST(CC, Diff2) {
	arma::mat a = { { 2, -1, 1 }, { -1, 2, 1 }, { 1, 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 2 } { -1 0 1 } { -2 -1 0 } }");
	CompatibleCartan cc;
	EXPECT_FALSE(cc(q, a));
}
}
