/**
 * cartan_iterator.h
 */
#pragma once
#ifndef REFL_CARTAN_ITERATOR_H__
#define REFL_CARTAN_ITERATOR_H__

#include <armadillo>

#include "qv/quiver_matrix.h"

namespace refl {
/**
 * Iterator over all possible signed quasi-Cartan matrices for a given quiver.
 */
class CartanIterator {
public:
	/**
	 * Make  new instance which provides quasi-Cartan matrices for the given
	 * quiver.
	 */
	CartanIterator(cluster::QuiverMatrix const& q);
	/**
	 * Check whether the iterator will return a valid matrix on the next call to
	 * next()
	 */
	bool
	has_next();
	/**
	 * Get the next Cartan matrix.
	 */
	arma::Mat<int>&
	next();
private:
	arma::Mat<int> _initial;
};
}
#endif

