/**
 * unique_cartan_filter.cc
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
#include "unique_cartan_filter.h"

namespace refl {
bool
UniqueCartanFilter::operator()(arma::Mat<int> const& m) {
	auto equiv= [&m,this](arma::Mat<int> const& other) {return _equiv_check(m, other);};
	bool not_found = std::find_if(_matrices.begin(), _matrices.end(), equiv) == _matrices.end();
	if(not_found) {
		_matrices.push_back(m);
	}
	return not_found;
}
}


