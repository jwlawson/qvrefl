/**
 * compatible_cartan.cc
 */
#include "compatible_cartan.h"

namespace refl {
bool
CompatibleCartan::operator()(cluster::QuiverMatrix const& lhs,
		arma::mat const& rhs) const {
	return false;
}
}

