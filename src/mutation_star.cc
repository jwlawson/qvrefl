/**
 * mutation_star.cc
 */
#include "mutation_star.h"

#include "util.h"

namespace refl {
MutationStar::MutationStar(cluster::QuiverMatrix const& q) 
	: _initial(q), 
		_qv_vector(q.num_cols(), cluster::QuiverMatrix(q.num_rows(), q.num_cols())),
		_mat_vector() {
	_mat_vector.reserve(q.num_cols());
	for(int32_t i = 0; i < q.num_cols(); ++i) {
		q.mutate(i, _qv_vector[i]);
		_mat_vector.push_back(util::to_arma(_qv_vector[i]));
	}
}

cluster::QuiverMatrix const&
MutationStar::qv(int32_t k) {
	return _qv_vector[k];
}
arma::Mat<int> const&
MutationStar::arma(int32_t k) {
	return _mat_vector[k];
}
}
