/**
 * unique_matrix_filter.h
 */
#pragma once
#ifndef REFL_UNIQUE_MATRIX_FILTER_H__
#define REFL_UNIQUE_MATRIX_FILTER_H__

#include <armadillo>
#include <unordered_set>

namespace refl {
class UniqueMatrixFilter {
	struct MatrixHash {
		std::size_t
		operator()(arma::Mat<int> const& m) const;
	};
	struct MatrixEqual {
		bool
		operator()(arma::Mat<int> const& lhs, arma::Mat<int> const& rhs) const;
	};
public:
	bool
	operator()(arma::Mat<int> const& m);
private:
	std::unordered_set<arma::Mat<int>, MatrixHash, MatrixEqual> _matrices;
};
}
#endif

