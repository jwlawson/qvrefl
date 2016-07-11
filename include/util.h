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
arma::mat
to_arma(cluster::QuiverMatrix const& q);
}
}
#endif
