/**
 * util.cc
 */
#include "util.h"

namespace refl {
namespace util {
arma::Mat<int>
to_arma(cluster::QuiverMatrix const& q) {
	/* arma matrices are stored in column first format, while QuiverMatrix is row
	 * first. The transpose is needed to fix this problem. */
	return arma::Mat<int>(q.data(), q.num_rows(), q.num_cols()).t();
}
}
}
