/**
 * semi_positive_filter.cc
 */
#include "semi_positive_filter.h"

namespace refl {
bool
SemiPositiveFilter::operator()(arma::Mat<int> const& a) {
	constexpr double tol = 1e-10;
	arma::mat d = arma::conv_to<arma::mat>::from(a);
	bool result = arma::eig_sym(_eigens, d);
	for(uint32_t i = 0; result && i < _eigens.n_elem; ++i) {
		result = _eigens(i) > -tol;
	}
	return result;
}
}
