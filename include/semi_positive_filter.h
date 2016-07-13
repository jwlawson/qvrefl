/**
 * semi_positive_filter.h
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
#ifndef REFL_SEMI_POSITIVE_FILTER_H__
#define REFL_SEMI_POSITIVE_FILTER_H__

#include <armadillo>

namespace refl {
class SemiPositiveFilter {
public:
	bool
	operator()(arma::Mat<int> const& a);
private:
	arma::vec _eigens;
	arma::mat _double_mat;
};
}
#endif
