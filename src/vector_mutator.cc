/**
 * vector_mutator.cc
 */
#include "vector_mutator.h"
namespace refl {
VectorMutator::VectorMutator(cluster::QuiverMatrix const& Q, arma::mat const& A)
	:	_quiver(Q),
		_cartan(A) {}
bool
VectorMutator::mutate(arma::mat const& vectors, uint32_t k, arma::mat& output) {
	for(int32_t i = 0; i < _quiver.num_rows(); ++i) {
		if(_quiver.get(i, k) > 0) {
			output.col(i) = vectors.col(i) - _cartan(i, k) * vectors.col(k);
		} else {
			output.col(i) = vectors.col(i);
		}
	}
	return true;
}
}

