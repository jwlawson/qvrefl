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
	for(uint_fast16_t i = 0; i < res.n_elem; ++i) {
		res(i) = std::abs(res(i));
	}
	for(uint_fast16_t i = 0; i < res.n_cols; ++i) {
		res(i, i) = 2;
	}
	return res;
}
/**
 * Make a bitset which has a 1 at each non-zero value in the off diagonal
 * matrix and a 0 at each zero.
 *
 * Then result.count() gives the number of non-zero entries, which is the number
 * of places the iterator needs to consider changing the signs of.
 *
 * @param a Initial quasi-Cartan matrix
 * @param vals The number of off diagonal upper triangle values in a
 */
boost::dynamic_bitset<>
zeros(arma::Mat<int> const& a, uint_fast16_t vals) {
	boost::dynamic_bitset<> result{ vals, 0 };
	uint_fast16_t row = 0;
	uint_fast16_t col = 1;
	for(uint_fast16_t i = 0; i < vals; ++i) {
		result[i] = a(row, col) != 0;
		if(++col >= a.n_cols) {
			col = ++row + 1;
		}
	}
	return result;
}
}
CartanIterator::CartanIterator(cluster::QuiverMatrix const& q) 
	: _number_vars( (q.num_rows() * (q.num_rows() - 1) ) / 2 ),
		_initial(convert(q)),
		_result(_initial.n_rows, _initial.n_cols),
		_zero_mask(zeros(_initial, _number_vars)),
		_non_zero_vars { _zero_mask.count() },
		_current_val {0},
		_max_val(std::pow(2, _non_zero_vars)) {}
bool
CartanIterator::has_next() {
	return _current_val < _max_val;
}
/*
 * Convert int _current_val into its binary representation, then insert this
 * into the non-zero values of _zero_mask. This ensures that we are never
 * pointlessly mutliplying 0 by -1 and creating many quasi-Cartans which are
 * actually the same.
 *
 * Use the bitset to multiply the vlues corresponding to 1s by -1 in the cartan
 * matrix.
 */
arma::Mat<int>&
CartanIterator::next(){
	boost::dynamic_bitset<> val_bits{_non_zero_vars, _current_val};
	boost::dynamic_bitset<> bits_with_zeros { _zero_mask };
	uint_fast16_t val_pos = 0;
	for(uint_fast16_t i = 0; i < bits_with_zeros.size(); ++i) {
		if(bits_with_zeros[i]) {
			bits_with_zeros[i] = val_bits[val_pos];
			++val_pos;
		}
	}

	_result = _initial;
	uint_fast16_t row = 0;
	uint_fast16_t col = 1;
	for(uint_fast16_t i = 0; i < _number_vars; ++i) {
		if(bits_with_zeros[i]) {
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

