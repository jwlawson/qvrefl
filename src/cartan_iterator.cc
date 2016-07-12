/**
 * cartan_iterator.cc
 */
#include "cartan_iterator.h"

#include "util.h"

#include <boost/dynamic_bitset.hpp>

namespace refl {
namespace {
arma::Mat<int>
convert(cluster::QuiverMatrix const& q) {
	arma::Mat<int> res(util::to_arma(q));
	for(uint32_t i = 0; i < res.n_elem; ++i) {
		res(i) = std::abs(res(i));
	}
	for(uint32_t i = 0; i < res.n_cols; ++i) {
		res(i, i) = 2;
	}
	return res;
}
}
CartanIterator::CartanIterator(cluster::QuiverMatrix const& q) 
	: _number_vars( (q.num_rows() * (q.num_rows() - 1) ) / 2 ),
		_initial(convert(q)),
		_result(_initial.n_rows, _initial.n_cols),
		_current_val {0},
		_max_val(std::pow(2, _number_vars)) {}
bool
CartanIterator::has_next() {
	return _current_val < _max_val;
}
arma::Mat<int>&
CartanIterator::next(){
	boost::dynamic_bitset<> bits{_number_vars, _current_val};
	_result = _initial;
	uint32_t row = 0;
	uint32_t col = 1;
	for(uint32_t i = 0; i < _number_vars; ++i) {
		if(bits[i]) {
			_result(row, col) = -1 * _result(row, col);
			_result(col, row) = -1 * _result(col, row);
		}
		if(++col >= _result.n_cols) {
			col = ++row + 1;
		}
	}
	++_current_val;
	return _result;
}
}

