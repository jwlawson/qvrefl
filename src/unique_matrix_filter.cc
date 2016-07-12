/**
 * unique_matrix_filter.cc
 */
#include "unique_matrix_filter.h"

namespace refl {
bool
UniqueMatrixFilter::operator()(arma::Mat<int> const& m) {
	return _matrices.insert(m).second;
}
std::size_t
UniqueMatrixFilter::MatrixHash::operator()(arma::Mat<int> const& m) const {
	std::size_t result = 17;
	uint64_t exp = 31;
	for(uint64_t i = 0; i < m.n_elem; ++i) {
		result *= exp;
		result += m(i);
	}
	return result;
}
bool
UniqueMatrixFilter::MatrixEqual::operator()(arma::Mat<int> const& lhs,
		arma::Mat<int> const& rhs) const {
	bool result = lhs.n_cols == rhs.n_cols && lhs.n_rows == rhs.n_rows;
	for(uint64_t i = 0; result && i < lhs.n_elem; ++i) {
		result = lhs(i) == rhs(i);
	}
	return result;
}
}

