/**
 * cartan.cc
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
#include "cartan_equiv.h"
#include "cartan_mutator.h"
#include "compatible_cartan.h"
#include "mutation_star.h"
#include "util.h"
#include "vector_mutator.h"

#include <fstream>
#include <iostream>
#include <unistd.h>


bool
check_compatible(cluster::QuiverMatrix const& q,
		arma::Mat<int> const& cartan, std::ostream& os = std::cout) {
	bool result = true;
	const refl::MutationStar star(q);
	const arma::mat initial_vecs(q.num_rows(), q.num_cols(), arma::fill::eye);
	arma::mat mutated_vecs(q.num_rows(), q.num_cols());
	arma::mat gram_matrix(q.num_rows(), q.num_cols());

	refl::CompatibleCartan compatible;
	if(!compatible(q, cartan) ) {
		result = false;
		os << "Initial matrix not compatible" << os.widen('\n');
	}
	refl::VectorMutator vmut(q, cartan);
	for(uint_fast16_t i = 0; i < cartan.n_cols; ++i) {
		vmut.mutate(initial_vecs, i, mutated_vecs);
		gram_matrix = refl::util::gram(mutated_vecs, cartan);
		if(!compatible(star.qv(i), gram_matrix) ) {
			result = false;
			os << "Mutation at " << i + 1 << " not compatible" << os.widen('\n');
		}
	}
	return result;
}
bool check_cartan(cluster::QuiverMatrix const& quiver, arma::Mat<int> const& cartan) {
	bool result = true;
	uint_fast16_t const nrows = cartan.n_rows;
	uint_fast16_t const ncols = cartan.n_cols;
	refl::CartanMutator cmut{quiver};
	refl::CartanEquiv cequiv;
	
	arma::Mat<int> mutated_cartan(nrows, ncols);
	cluster::QuiverMatrix mutated_quiver(nrows, ncols);

	for (uint_fast16_t k = 0; k < cartan.n_rows; ++k) {
		cmut(cartan, k, mutated_cartan);
		quiver.mutate(k, mutated_quiver);

		if(!check_compatible(mutated_quiver, mutated_cartan) ) {
			result = false;
			std::cout << k + 1 << ": " << "not compatible" << std::cout.widen('\n')
				<< "Mutated quiver:\n" << mutated_quiver << std::cout.widen('\n');
			mutated_cartan.print("Mutated cartan:");
		}
	}
	return result;
}
void run_on_input(std::istream& is) {
	std::string line;
	arma::Mat<int> c;
	size_t count = 0;
	while(is.good() ) {
		cluster::QuiverMatrix *q = nullptr;
		if(std::getline(is, line) && line.length() > 4 ) {
			q = new cluster::QuiverMatrix(line);
		} 
		if( q && std::getline(is, line) && line.length() > 4 ) {
			cluster::QuiverMatrix cq {line};
			std::cout << ++count << ": " << *q << std::cout.widen('\n');
			c = refl::util::to_arma(cq);
			check_compatible(*q, c);
		}
		if( q && !c.empty() ) {
			check_cartan( *q, c );
		}
		if( q ) { delete q; }
	}
}

void usage() {
	std::cout << "cartan -i in_file" << '\n';
}
int main(int argc, char* argv[]) {
	std::string in_filename;
	int c;
	while ((c = getopt (argc, argv, "i:")) != -1) {
    switch (c) {
			case 'i':
				in_filename = optarg;
				break;
			case 'h':
      case '?':
      default:
				usage();
				return 1;
		}
	}
	if( in_filename.empty() ) {
		usage();
		return 1;
	}
	std::ifstream in_file;
	in_file.open(in_filename);
	if(! in_file.is_open()) {
		std::cerr << "Cannot open file " << in_filename << std::endl;
		return 2;
	}
	run_on_input(in_file);
	return 0;
}

