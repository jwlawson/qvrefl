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
	for(int_fast16_t i = 0; i < q.num_cols(); ++i) {
		q.mutate(i, _qv_vector[i]);
		_mat_vector.push_back(util::to_arma(_qv_vector[i]));
	}
}

cluster::QuiverMatrix const&
MutationStar::qv(uint_fast16_t k) const {
	return _qv_vector[k];
}
arma::Mat<int> const&
MutationStar::arma(uint_fast16_t k) const {
	return _mat_vector[k];
}
}
