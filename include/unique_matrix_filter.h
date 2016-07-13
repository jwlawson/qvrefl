/**
 * unique_matrix_filter.h
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

