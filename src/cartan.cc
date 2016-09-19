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
#include "util.h"
#include "vector_mutator.h"

#include <fstream>
#include <iostream>
#include <unistd.h>


void check_cartan(cluster::QuiverMatrix const& quiver, arma::Mat<int> const& cartan) {
	uint_fast16_t const nrows = cartan.n_rows;
	uint_fast16_t const ncols = cartan.n_cols;
	arma::Mat<int> const vecs(nrows, ncols, arma::fill::eye);
	refl::VectorMutator vmut{quiver, cartan};
	refl::CartanMutator cmut{quiver};
	refl::CartanEquiv cequiv;
	
	arma::Mat<int> mutated_vecs(nrows, ncols);
	arma::Mat<int> mutated_cartan(nrows, ncols);
	arma::Mat<int> gram(nrows, ncols);

	std::cout << quiver << std::cout.widen('\n');

	for (uint_fast16_t k = 0; k < cartan.n_rows; ++k) {
		vmut.mutate(vecs, k, mutated_vecs);
		gram = refl::util::gram( mutated_vecs, cartan );
		cmut(cartan, k, mutated_cartan);

		if(!cequiv(mutated_cartan, gram) ) {
			std::cout << k << ": " << "not compatible" << std::cout.widen('\n');
			mutated_cartan.print("Mutated cartan:");
			gram.print("Gram of mutated vectors:");
		}
	}
}
void run_on_input(std::istream& is) {
	std::string line;
	arma::Mat<int> c;
	while(is.good() ) {
		cluster::QuiverMatrix *q = nullptr;
		if(std::getline(is, line) && line.length() > 4 ) {
			q = new cluster::QuiverMatrix(line);
		} 
		if( std::getline(is, line) && line.length() > 4 ) {
			cluster::QuiverMatrix cq {line};
			c = refl::util::to_arma(cq);
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

