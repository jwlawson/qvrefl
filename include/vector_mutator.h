/**
 * vector_mutator.h
 */
#pragma once
#ifndef REFL_VECTOR_MUTATOR_H__
#define REFL_VECTOR_MUTATOR_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
class VectorMutator {
public:
	/**
	 * Construct a mutator with quasi-Cartan matrix A for quiver Q.
	 */
	VectorMutator(cluster::QuiverMatrix const& Q, arma::mat const& A);
	/**
	 * Mutate the vectors in the direction k, using the set quasi-Cartan matrix
	 * and quiver. The resulting vectors are computed in the provided output matrix.
	 */
	bool
	mutate(arma::mat const& vectors, uint32_t k, arma::mat& output);
};
}
#endif

