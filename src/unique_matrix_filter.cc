/**
 * unique_matrix_filter.cc
 */
#include "unique_matrix_filter.h"

namespace refl {
bool
UniqueMatrixFilter::operator()(arma::Mat<int> const& m) {
	return false;
}
}

