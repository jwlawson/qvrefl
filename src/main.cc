/**
 * main.cc
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
#include <iostream>
#include <string>
#include <unistd.h>

#include <boost/dynamic_bitset.hpp>

#include "filtered_iterator.h"
#include "cartan_iterator.h"
#include "compatible_cartan.h"
#include "mutation_star.h"
#include "semi_positive_filter.h"
#include "vector_mutator.h"
#include "unique_matrix_filter.h"
#include "util.h"

namespace refl {
namespace {
using MatrixVec = std::vector<arma::Mat<int>>;
using UniqueCartanIter = FilteredIterator<CartanIterator, arma::Mat<int>, UniqueMatrixFilter>;
using SemiPosIter = FilteredIterator<UniqueCartanIter, arma::Mat<int>, SemiPositiveFilter>;

SemiPosIter
get_cartan_iterator(cluster::QuiverMatrix const& q) {
	CartanIterator cartan(q);
	UniqueCartanIter unique(std::move(cartan));
	return SemiPosIter(std::move(unique));
}
MatrixVec
find_serving_cartans(cluster::QuiverMatrix const& q) {
	MatrixVec result(q.num_cols(), arma::Mat<int>(q.num_rows(), q.num_cols()));
	boost::dynamic_bitset<> found{ q.num_cols(), 0 };
	const MutationStar star(q);
	const arma::mat initial_vecs(q.num_rows(), q.num_cols(), arma::fill::eye);
	arma::mat mutated_vecs(q.num_rows(), q.num_cols());
	arma::mat gram_matrix(q.num_rows(), q.num_cols());

	CompatibleCartan compatible;

	SemiPosIter pos_iter = get_cartan_iterator(q);

	while(pos_iter.has_next() && !found.all()) {
		arma::Mat<int> const& AQ = pos_iter.next();
		VectorMutator vmut(q, AQ);
		for(uint_fast16_t i = 0; i < AQ.n_cols; ++i) {
			vmut.mutate(initial_vecs, i, mutated_vecs);
			gram_matrix = util::gram(mutated_vecs, AQ);
			if(!found[i] && compatible(star.qv(i), gram_matrix)) {
				result[i] = AQ;
				found[i] = true;
			}
		}
	}

	return result;
}
}
}
void
usage() {
	std::cout << "qvrefl -m matrix" << std::cout.widen('\n');
	std::cout << "  -m Specify matrix to find quasi-Cartan companions of" << std::cout.widen('\n');
	std::cout.flush();
}
int
main(int argc, char* argv[]) {
	std::string matrix;
	int c;
	while ((c = getopt (argc, argv, "m:h")) != -1) {
    switch (c) {
      case 'm':
				matrix = optarg;
        break;
			case 'h':
      case '?':
      default:
				usage();
				return 1;
		}
	}
	if(matrix.length() < 4) {
		usage();
		return 2;
	}
	cluster::QuiverMatrix q(matrix);
	auto res = refl::find_serving_cartans(q);
	for(uint_fast16_t i = 0; i < res.size(); ++i) {
		if(res[i].at(0, 0) != 2) {
			std::cout << i << " is not served by any semi-positive quasi-Cartan" << std::cout.widen('\n');
		} else {
			std::string label;
			label.append(std::to_string(i)).append(" served by");
			res[i].print(label);
		}
	}
	std::cout.flush();
	return 0;
}

