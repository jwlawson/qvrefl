#include "fully_compatible_check.h"

#include "cartan_mutator.h"
#include "compatible_cartan.h"

namespace refl {
bool FullyCompatibleCheck::operator()(cluster::QuiverMatrix const& q,
                                      MutationStar const& star,
                                      arma::Mat<int> const& AQ) {
  static arma::Mat<int> gram_matrix;
  bool result = true;
  gram_matrix.set_size(q.num_rows(), q.num_cols());

  CompatibleCartan compatible;
  result = compatible(q, AQ);

  CartanMutator cmut(q);
  int_fast16_t ncols = AQ.n_cols;
  for (int_fast16_t i = 0; result && i < ncols; ++i) {
    cmut(AQ, i, gram_matrix);
    result = compatible(star.qv(i), gram_matrix);
  }
  return result;
}
bool FullyCompatibleCheck::operator()(cluster::QuiverMatrix const& q,
                                      arma::Mat<int> const& AQ) {
  const MutationStar star(q);
  return operator()(q, star, AQ);
}
}
