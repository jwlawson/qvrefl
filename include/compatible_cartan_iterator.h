#pragma once
#ifndef _REFL_COMPATIBLE_CARTAN_ITER_H__
#define _REFL_COMPATIBLE_CARTAN_ITER_H__

#include "cartan_iterator.h"
#include "filtered_iterator.h"
#include "fully_compatible_check.h"
#include "semi_positive_filter.h"
#include "unique_cartan_filter.h"
#include "unique_matrix_filter.h"

#include "qv/quiver_matrix.h"

#include <armadillo>

namespace refl {
namespace compatible_cartan_iter {
using UniqueCartanIter = FilteredIterator<CartanIterator, arma::Mat<int>, UniqueMatrixFilter>;
using NonEquivCartanIter = FilteredIterator<UniqueCartanIter, arma::Mat<int>, UniqueCartanFilter>;
using CompatibleIter = FilteredIterator<NonEquivCartanIter, arma::Mat<int>, FixedQuiverFullyCompatibleCheck>;
using SemiPosIter = FilteredIterator<CompatibleIter, arma::Mat<int>, SemiPositiveFilter>;
SemiPosIter Get(cluster::QuiverMatrix const& q) {
	CartanIterator ci(q);
	UniqueCartanIter un(std::move(ci));
	NonEquivCartanIter nequ(std::move(un));
	FixedQuiverFullyCompatibleCheck chk(q);
	CompatibleIter comp(std::move(nequ), std::move(chk));
	return SemiPosIter(std::move(comp));
}
}
/**
 * An iterator which provides semipositive compatible quasi-cartan companions
 * for a given quiver.
 */
struct CompatibleCartanIterator : compatible_cartan_iter::SemiPosIter {
	CompatibleCartanIterator(cluster::QuiverMatrix const& q) :
		compatible_cartan_iter::SemiPosIter(compatible_cartan_iter::Get(q)) {}
	using compatible_cartan_iter::SemiPosIter::has_next;
	using compatible_cartan_iter::SemiPosIter::next;
};
}
#endif
