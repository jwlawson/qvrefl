/**
 * util.cc
 */
#include "util.h"

namespace refl {
namespace util {
arma::Mat<int>
to_arma(cluster::QuiverMatrix const& q) {
	return arma::Mat<int>(q.data(), q.num_rows(), q.num_cols());
}
}
}
