/**
 * util.cc
 */
#include "util.h"

namespace refl {
namespace util {
bool
equal(arma::mat const& a, arma::mat const& b) {
	constexpr double tol = 1e-10;
	bool result = a.n_cols == b.n_cols && a.n_rows == b.n_rows;
	for(auto a_it = a.begin(), a_end = a.end(), b_it = b.begin(), b_end = b.end();
			a_it != a_end && b_it != b_end;
			++a_it, ++b_it) {
		result = result && std::abs(*a_it - *b_it) < tol;
	}
	return result;
}
arma::mat
to_arma(cluster::QuiverMatrix const& q) {
	return arma::mat();
}
arma::mat
gram(arma::mat const& vectors, arma::mat const& product){
	return arma::mat();
}
}
}
