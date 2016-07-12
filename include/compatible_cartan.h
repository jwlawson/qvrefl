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
	template<class mat>
	bool
	operator()(cluster::QuiverMatrix const& lhs, mat const& rhs) const;
};
template<class mat>
bool
CompatibleCartan::operator()(cluster::QuiverMatrix const& lhs,
		mat const& rhs) const {
	bool result = static_cast<uint32_t>(lhs.num_cols()) == rhs.n_cols
		&& static_cast<uint32_t>(lhs.num_rows()) == rhs.n_rows;
	/* As the two types of matrix are arranged differently in memory, we cannot
	 * use quick comparisons, but have to go through each value separately. */
	for(uint32_t row = 0; result && row < rhs.n_rows; ++row) {
		for(uint32_t col = 0; result && col < rhs.n_cols; ++col) {
			if(row == col) continue;
			result = std::abs(lhs.get(row,col)) == std::abs(rhs(row, col));
		}
	}
	return result;
}
}
#endif

