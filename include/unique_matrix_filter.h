/**
 * unique_matrix_filter.h
 */
#pragma once
#ifndef REFL_UNIQUE_MATRIX_FILTER_H__
#define REFL_UNIQUE_MATRIX_FILTER_H__

#include <armadillo>

namespace refl {
class UniqueMatrixFilter {
public:
	bool
	operator()(arma::Mat<int> const& m);
};
}
#endif

