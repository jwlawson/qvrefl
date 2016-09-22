/**
 * qvrefl.cc
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

#include "qv/mutation_class_loader.h"
#include "qv/equiv_mutation_class_loader.h"
#include "qv/stream_iterator.h"

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
void
output_serving_cartans(cluster::QuiverMatrix const& q, std::ostream& os = std::cout) {
	auto res = refl::find_serving_cartans(q);
	for(uint_fast16_t i = 0; i < res.size(); ++i) {
		if(res[i].at(0, 0) != 2) {
			os << i << " is not served by any semi-positive quasi-Cartan" << os.widen('\n');
		} else {
			std::string label;
			label.append(std::to_string(i)).append(" served by");
			res[i].print(os, label);
		}
	}
}
std::pair<bool, arma::Mat<int> >
has_same_serving_cartan(cluster::QuiverMatrix const& q) {
	bool result = false;
	arma::Mat<int> cartan;
	const MutationStar star(q);
	const arma::mat initial_vecs(q.num_rows(), q.num_cols(), arma::fill::eye);
	arma::mat mutated_vecs(q.num_rows(), q.num_cols());
	arma::mat gram_matrix(q.num_rows(), q.num_cols());

	CompatibleCartan compatible;

	SemiPosIter pos_iter = get_cartan_iterator(q);

	while(pos_iter.has_next() && !result) {
		bool is_comp = true;
		arma::Mat<int> const& AQ = pos_iter.next();
		VectorMutator vmut(q, AQ);
		for(uint_fast16_t i = 0; is_comp && i < AQ.n_cols; ++i) {
			vmut.mutate(initial_vecs, i, mutated_vecs);
			gram_matrix = util::gram(mutated_vecs, AQ);
			is_comp = compatible(star.qv(i), gram_matrix);
		}
		result = is_comp;
		if(result) {
			cartan = AQ;
		}
	}
	return {result, cartan};
}
void
check_compatible(cluster::QuiverMatrix const& q,
		cluster::QuiverMatrix const& cartan, std::ostream& os = std::cout) {
	arma::Mat<int> AQ = util::to_arma(cartan);
	const MutationStar star(q);
	const arma::mat initial_vecs(q.num_rows(), q.num_cols(), arma::fill::eye);
	arma::mat mutated_vecs(q.num_rows(), q.num_cols());
	arma::mat gram_matrix(q.num_rows(), q.num_cols());

	CompatibleCartan compatible;
	if(! compatible(q, AQ) ) {
		os << "Initial matrix not compatible" << os.widen('\n');
	}
	VectorMutator vmut(q, AQ);
	for(uint_fast16_t i = 0; i < AQ.n_cols; ++i) {
		vmut.mutate(initial_vecs, i, mutated_vecs);
		gram_matrix = util::gram(mutated_vecs, AQ);
		if(!compatible(star.qv(i), gram_matrix) ) {
			os << "Mutation at " << i << " not compatible" << os.widen('\n');
		}
	}
}
}
}
enum Func {
	ListCompatible,
	SameCompatible,
	CheckSingle
};
void
usage() {
	std::cout << "qvrefl -cls [-m matrix] [-i in_file] [-a cartan]" << std::cout.widen('\n');
	std::cout << "  -m Specify matrix to find quasi-Cartan companions of" << std::cout.widen('\n');
	std::cout << "  -i Specify input file of matrices to read" << std::cout.widen('\n');
	std::cout << "  -c Check the whole mutation class of the matrix" << std::cout.widen('\n');
	std::cout << "  -l List the first quasi-Cartan serving each vertex" << std::cout.widen('\n');
	std::cout << "  -s Check if there is a quasi-Cartan serving all vertices" << std::cout.widen('\n');
	std::cout << "  -a Specify a cartan matrix to check whether it serves every vertex" << std::cout.widen('\n');
	std::cout.flush();
}
int
main(int argc, char* argv[]) {
	std::string matrix;
	std::string cartan;
	Func function = Func::ListCompatible;
	std::string input;
	bool mut_class = false;
	int c;
	while ((c = getopt (argc, argv, "lsm:i:hca:")) != -1) {
    switch (c) {
      case 'm':
				matrix = optarg;
        break;
			case 'a':
				cartan = optarg;
				function = Func::CheckSingle;
				break;
			case 'l':
				function = Func::ListCompatible;
				break;
			case 's':
				function = Func::SameCompatible;
				break;
			case 'i':
				input = optarg;
				break;
			case 'c':
				mut_class = true;
				break;
			case 'h':
      case '?':
      default:
				usage();
				return 1;
		}
	}
	if(matrix.length() < 4 && input.empty()) {
		usage();
		return 2;
	}
	if(function == Func::ListCompatible) {
		if(matrix.length() > 0) {
			if(mut_class) {
				cluster::EquivQuiverMatrix q(matrix);
				cluster::EquivMutationClassLoader cl(q);
				while(cl.has_next()) {
					auto quiver = cl.next_ptr();
					std::cout << *quiver << ":" << std::cout.widen('\n');
					refl::output_serving_cartans(*quiver);
				}
			} else {
				cluster::QuiverMatrix q(matrix);
				refl::output_serving_cartans(q);
			}
		} else {
			std::ifstream inf;
			inf.open(input);
			if(!inf.is_open()) {
				std::cerr << "Could not open file " << input << std::endl;
				return 1;
			}
			cluster::StreamIterator<cluster::QuiverMatrix> iter(inf);
			while(iter.has_next()) {
				auto quiver = iter.next();
				std::cout << *quiver << std::cout.widen('\n');
				refl::output_serving_cartans(*quiver);
			}
		}
	} else if (function == Func::SameCompatible) {
		if(matrix.length() > 0) {
			if(mut_class) {
				cluster::EquivQuiverMatrix q(matrix);
				cluster::EquivMutationClassLoader cl(q);
				while(cl.has_next()) {
					auto quiver = cl.next_ptr();
					bool result = refl::has_same_serving_cartan(*quiver).first;
					std::cout << (result ? "True: " : "False: ") << *quiver << std::cout.widen('\n');
				}
			} else {
				cluster::QuiverMatrix q(matrix);
				auto result = refl::has_same_serving_cartan(q);
				std::cout << (result.first ? result.second : "False") << std::cout.widen('\n');
			}
		} else {
			std::ifstream inf;
			inf.open(input);
			if(!inf.is_open()) {
				std::cerr << "Could not open file " << input << std::endl;
				return 1;
			}
			cluster::StreamIterator<cluster::QuiverMatrix> iter(inf);
			while(iter.has_next()) {
				auto quiver = iter.next();
				bool result = refl::has_same_serving_cartan(*quiver).first;
				std::cout << (result ? "True: " : "False: ") << *quiver << std::cout.widen('\n');
			}
		}
	} else if (function == Func::CheckSingle ) {
		if(matrix.empty() || cartan.empty() ) {
			usage();
			return 4;
		}
		refl::check_compatible( cluster::QuiverMatrix{matrix}, cluster::QuiverMatrix{cartan} );
	} else {
		usage();
	}
	std::cout.flush();
	return 0;
}

