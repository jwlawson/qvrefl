#include "cartan_exchange_graph.h"

#include "cartan_mutator.h"
#include "fully_compatible_check.h"
#include "permuted_cartan_equiv.h"

namespace refl {
namespace cartan_exchange {
void
Mutator::operator()(CartanQuiver const* const initial, size_t k,
                    CartanQuiver& output) const {
  initial->quiver.mutate(k, output.quiver);
  CartanMutator cmut(initial->quiver);
  cmut(initial->cartan, k, output.cartan);
  FixedQuiverFullyCompatibleCheck comp(output.quiver);
  output.fully_compatible = comp(output.cartan);
}
bool
Equiv::operator()(CartanQuiver const* const lhs,
                  CartanQuiver const* const rhs) const {
  auto check = cluster::EquivalenceChecker::Get(lhs->quiver.num_rows());
  bool result = check->are_equivalent(lhs->quiver, rhs->quiver);
  if (result) {
    auto permutation = check->last_row_map();
    PermutedCartanEquiv equiv;
    result = equiv(lhs->cartan, rhs->cartan, permutation);
    if (!result) {
      // valid_row_maps assumes the matrices are equivalent, so we do need the
      // equivalent check above before we can check all permutations.
      auto all_perms = check->valid_row_maps(lhs->quiver, rhs->quiver);
      for (auto permutation : *all_perms) {
        result = equiv(lhs->cartan, rhs->cartan, permutation);
        if (result) {
          break;
        }
      }
    }
  }
  return result;
}
bool
Equal::operator()(CartanQuiver const* const lhs,
                  CartanQuiver const* const rhs) const {
  return cluster::IntMatrix::are_equal(lhs->quiver, rhs->quiver);
}
size_t
Hash::operator()(CartanQuiver const* const quiver) const {
  return quiver->quiver.hash();
}
CartanQuiver*
NewInstance::operator()(size_t const size) const {
  return new CartanQuiver{cluster::EquivQuiverMatrix(size, size),
                          arma::Mat<int>(size, size), false};
}
bool
DontMutateNonCompatible::operator()(CartanQuiver const* const mptr,
                                    int /*vertex*/) const {
  return mptr->fully_compatible;
}
}  // namespace cartan_exchange
}  // namespace refl
