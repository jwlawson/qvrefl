/**
 * util.h
 * Copyright 2016 John Lawson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#ifndef REFL_UTIL_H__
#define REFL_UTIL_H__

#include <armadillo>
#include "qv/quiver_matrix.h"

namespace refl {
namespace util {
/**
 * Check if two arma::mats are equal.
 */
template<class mat>
bool
equal(mat const& lhs, mat const& rhs);
/**
 * Convert a cluster QuiverMatrix to an armadillo matrix.
 */
arma::Mat<int>
to_arma(cluster::QuiverMatrix const& q);
/**
 * Construct the Gram matrix of the specified vectors using the provided
 * symmetric matrix as the inner product.
 */
template<class mat1, class mat2>
mat1
gram(mat1 const& vectors, mat2 const& prod);
}

template<class mat>
bool
util::equal(mat const& a, mat const& b) {
	constexpr double tol = 1e-10;
	bool result = a.n_cols == b.n_cols && a.n_rows == b.n_rows;
	for(auto a_it = a.begin(), a_end = a.end(), b_it = b.begin(), b_end = b.end();
			a_it != a_end && b_it != b_end;
			++a_it, ++b_it) {
		result = result && std::abs(*a_it - *b_it) < tol;
	}
	return result;
}
template<class mat1, class mat2>
mat1
util::gram(mat1 const& vectors, mat2 const& product){
	return vectors.t() * product * vectors;
}
}
#endif
