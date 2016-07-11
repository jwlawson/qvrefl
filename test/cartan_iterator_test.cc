/**
 * cartan_iterator_test.cc
 */
#include "cartan_iterator.h"

#include "gtest/gtest.h"

#include "util.h"

namespace refl {
TEST(CartanIt, 2Dim) {
	cluster::QuiverMatrix q("{ { 0 1 } { -1 0 } }");
	CartanIterator c(q);
	ASSERT_TRUE(c.has_next());
	arma::mat& a = c.next();
	arma::mat exp = { { 2, 1 }, { 1, 2 } };
	EXPECT_TRUE(util::equal(a, exp));
	ASSERT_TRUE(c.has_next());
	a = c.next();
	arma::mat exp2 = { { 2, -1 }, { -1, 2 } };
	EXPECT_TRUE(util::equal(a, exp2));
	ASSERT_FALSE(c.has_next());
}
TEST(CartanIt, Count4) {
	cluster::QuiverMatrix q("{ { 0 1 0 0 } { -1 0 1 0 } { 0 -1 0 1 } { 0 0 -1 0 } }");
	CartanIterator c(q);
	ASSERT_TRUE(c.has_next());
	int count = 0;
	while(c.has_next()) {
		c.next();
		++count;
	}
	EXPECT_EQ(16, count);
}
}
