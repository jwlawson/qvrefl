/**
 * vector_mutator.cc
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
#include "vector_mutator.h"
namespace refl {
VectorMutator::VectorMutator(cluster::QuiverMatrix const& Q, arma::Mat<int> const& A)
	:	_quiver(Q),
		_cartan(A) {}
bool
VectorMutator::mutate(arma::mat const& vectors, uint_fast16_t k, arma::mat& output) {
	for(int_fast16_t i = 0; i < _quiver.num_rows(); ++i) {
		if(_quiver.get(i, k) > 0) {
			output.col(i) = vectors.col(i) - _cartan(i, k) * vectors.col(k);
		} else {
			output.col(i) = vectors.col(i);
		}
	}
	return true;
}
}

