/*
 * filtered_iterator.h
 * Copyright 2016 John Lawson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
inline
FilteredIterator<It, Output, Filter, positive>::FilteredIterator(It && it)
	: _it(std::move(it)),
		_result(),
		_next(),
		_filter(),
		_got_next{true} {
	get_next();
}
template <class It, class Output, class Filter, bool positive>
inline
bool
FilteredIterator<It, Output, Filter, positive>::has_next() {
	return _got_next;
}
template <class It, class Output, class Filter, bool positive>
inline
const Output & 
FilteredIterator<It, Output, Filter, positive>::next() {
	_result.swap(_next);
	get_next();
	return _result;
}
template <class It, class Output, class Filter, bool positive>
inline
void 
FilteredIterator<It, Output, Filter, positive>::get_next() {
	if(!_it.has_next()) {
		_got_next = false;
		return;
	}
	do {
		_next = _it.next();
	} while(_filter(_next) != positive && _it.has_next());

	_got_next = _it.has_next() || _filter(_next) == positive;
}
}
#endif

