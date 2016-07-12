/**
 * semi_positive_filter.h
 */
#pragma once
#ifndef REFL_SEMI_POSITIVE_FILTER_H__
#define REFL_SEMI_POSITIVE_FILTER_H__

#include <armadillo>

namespace refl {
class SemiPositiveFilter {
public:
	bool
	operator()(arma::Mat<int> const& a);
private:
	arma::vec _eigens;
};
}
#endif
