/*
 * filtered_iterator.h
 */
#pragma once
#ifndef REFL_FILTERED_ITERATOR_H__
#define REFL_FILTERED_ITERATOR_H__

#include <utility>

namespace refl {
/**
 * Filters the output of the given iterator.
 * Only those outputs where Filter(output) == positive are outputted.
 */
template <class It, class Output, class Filter, bool positive = true>
class FilteredIterator {
public:
	/** Create an iterator which filters the provided iterator. */
	FilteredIterator(It && it);
	/**
	 * Check whether the iterator will return a valid output on the next call of
	 * next()
	 */
	bool
	has_next();
	/**
	 * Get the next filtered output of the iterator.
	 */
	Output const&
	next();
private:
	It _it;
	Output _result;
	Output _next;
	Filter _filter;
	bool _got_next;

	/** Compute the next value to return. */
	void
	get_next();
};
template <class It, class Output, class Filter, bool positive>
FilteredIterator<It, Output, Filter, positive>::FilteredIterator(It && it)
	: _it(std::move(it)),
		_result(),
		_next(),
		_filter(),
		_got_next{true} {
	get_next();
}
template <class It, class Output, class Filter, bool positive>
bool
FilteredIterator<It, Output, Filter, positive>::has_next() {
	return _got_next;
}
template <class It, class Output, class Filter, bool positive>
const Output & 
FilteredIterator<It, Output, Filter, positive>::next() {
	_result.swap(_next);
	get_next();
	return _result;
}
template <class It, class Output, class Filter, bool positive>
void 
FilteredIterator<It, Output, Filter, positive>::get_next() {
	if(!_it.has_next()) {
		return;
	}
	do {
		_next = _it.next();
	} while(_filter(_next) != positive && _it.has_next());

	_got_next = _it.has_next() || _filter(_next) == positive;
}
}
#endif

