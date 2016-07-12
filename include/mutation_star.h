/**
 * mutation_star.h
 */
#pragma once
#ifndef REFL_MUTATION_STAR_H__
#define REFL_MUTATION_STAR_H__

#include "qv/quiver_matrix.h"

#include <armadillo>

namespace refl {
/**
 * Container of all quivers obtainable from the provided quiver by a single
 * mutation.
 *
 * These quivers can be accessed as either QuiverMatrix or arma::mat objects.
 */
class MutationStar {
public:
	/**
	 * Construct a new instance with the provided initial quiver.
	 */
	MutationStar(cluster::QuiverMatrix const& q);
	/**
	 * Get the mutation in the k-th direction of the initial quiver.
	 */
	cluster::QuiverMatrix const&
	qv(int32_t k);
	arma::Mat<int> const&
	arma(int32_t k);
private:
	cluster::QuiverMatrix const& _initial;
	std::vector<cluster::QuiverMatrix> _qv_vector;
	std::vector<arma::Mat<int>> _mat_vector;
};
}
#endif

