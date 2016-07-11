/**
 * util.h
 */
#pragma once
#ifndef REFL_UTIL_H__
#define REFL_UTIL_H__

#include <armadillo>
#include "qv/quiver_matrix.h"

namespace refl {
namespace util {
/**
 * Convert a cluster QuiverMatrix to an armadillo matrix.
 */
arma::mat
to_arma(cluster::QuiverMatrix const& q);
/**
 * Construct the Gram matrix of the specified vectors using the provided
 * symmetric matrix as the inner product.
 */
arma::mat
gram(arma::mat const& vectors, arma::mat const& prod);
}
}
#endif
