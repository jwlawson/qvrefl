/**
 * semi_positive_filter.cc
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
#include "semi_positive_filter.h"

namespace refl {
bool
SemiPositiveFilter::operator()(arma::Mat<int> const& a) {
	constexpr float tol = 1e-5;
	_mat = arma::conv_to<arma::Mat<float>>::from(a);
	bool result = arma::eig_sym(_eigens, _mat);
	for(uint_fast16_t i = 0; result && i < _eigens.n_elem; ++i) {
		result = _eigens(i) > -tol;
	}
	return result;
}
}
