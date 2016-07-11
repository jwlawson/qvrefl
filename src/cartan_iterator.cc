/**
 * cartan_iterator.cc
 */
#include "cartan_iterator.h"

#include "util.h"

namespace refl {
CartanIterator::CartanIterator(cluster::QuiverMatrix const& q) 
	: _initial(util::to_arma(q)) {}
bool
CartanIterator::has_next() {
	return false;
}
arma::Mat<int>&
CartanIterator::next(){
	return _initial;
}
}

