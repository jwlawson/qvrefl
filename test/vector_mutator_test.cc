/**
 * vector_mutator_test.cc
 */
#include "vector_mutator.h"

#include "gtest/gtest.h"
#include "qv/quiver_matrix.h"

#include <armadillo>

namespace refl {
TEST(VecMut, A2) {
	arma::mat a = { { 2, 1 }, { 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 } { -1 0 } }");
	VectorMutator vm(q, a);
	arma::mat vecs = { { 0, 1 }, { 1, 0 } };
	arma::mat res(2, 2);
	vm.mutate(vecs, 0, res);
	arma::mat exp = { { 0, 1 }, { 1, 0 } };
	EXPECT_TRUE(arma::all(arma::all(exp == res)));

	arma::mat res2(2, 2);
	vm.mutate(vecs, 1, res2);
	arma::mat exp2 = { { -1, 1 }, { 1, 0 } };
	EXPECT_TRUE(arma::all(arma::all(exp2 == res2)));
}
TEST(VecMut, A3) {
	arma::mat a = { { 2, 1, -1 }, { 1, 2, 1 }, { -1, 1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 0 } { -1 0 1 } { 0 -1 0 } }");
	VectorMutator vm(q, a);
	arma::mat vecs = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	arma::mat res(3, 3);
	vm.mutate(vecs, 0, res);
	arma::mat exp = vecs;
	EXPECT_TRUE(arma::all(arma::all(exp == res)));

	arma::mat res2(3, 3);
	vm.mutate(vecs, 1, res2);
	arma::mat exp2 = { { 1, 0, 0 }, { -1, 1, 0 }, { 0, 0, 1 } };
	EXPECT_TRUE(arma::all(arma::all(exp2 == res2)));

	arma::mat res3(3, 3);
	vm.mutate(vecs, 2, res3);
	arma::mat exp3 = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 1 } };
	EXPECT_TRUE(arma::all(arma::all(exp2 == res2)));
}
TEST(VecMut, A3signs) {
	arma::mat a = { { 2, -1, -1 }, { -1, 2, -1 }, { -1, -1, 2 } };
	cluster::QuiverMatrix q("{ { 0 1 0 } { -1 0 1 } { 0 -1 0 } }");
	VectorMutator vm(q, a);
	arma::mat vecs = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	arma::mat res(3, 3);
	vm.mutate(vecs, 0, res);
	arma::mat exp = vecs;
	EXPECT_TRUE(arma::all(arma::all(exp == res)));

	arma::mat res2(3, 3);
	vm.mutate(vecs, 1, res2);
	arma::mat exp2 = { { 1, 0, 0 }, { 1, 1, 0 }, { 0, 0, 1 } };
	EXPECT_TRUE(arma::all(arma::all(exp2 == res2)));

	arma::mat res3(3, 3);
	vm.mutate(vecs, 2, res3);
	arma::mat exp3 = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 1, 1 } };
	EXPECT_TRUE(arma::all(arma::all(exp2 == res2)));
}
}
