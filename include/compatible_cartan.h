/**
 * compatible_cartan.h
 */
#pragma once
#ifndef REFL_COMPATIBLE_CARTAN_H__
#define REFL_COMPATIBLE_CARTAN_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
class CompatibleCartan {
public:
	bool
	operator()(cluster::QuiverMatrix const& lhs, arma::mat const& rhs) const;
};
}
#endif

