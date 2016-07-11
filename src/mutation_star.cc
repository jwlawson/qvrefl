/**
 * mutation_star.cc
 */
#include "mutation_star.h"

namespace refl {
MutationStar::MutationStar(cluster::QuiverMatrix const& q) 
	: _initial(q) {}
cluster::QuiverMatrix const&
MutationStar::qv(int32_t k) {
	return _initial;
}
arma::mat const&
MutationStar::arma(int32_t k) {
	return _fake;
}
}
