/**
 * compatible_cartan.h
 */
#pragma once
#ifndef REFL_COMPATIBLE_CARTAN_H__
#define REFL_COMPATIBLE_CARTAN_H__

#include <armadillo>

namespace refl {
class CompatibleCartan {
public:
	bool
	operator()(arma::mat const& lhs, arma::mat const& rhs) const;
};
}
#endif

