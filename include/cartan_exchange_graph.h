#pragma once
#ifndef REFL_CARTAN_EXCHANGE_GRAPH_H__
#define REFL_CARTAN_EXCHANGE_GRAPH_H__

#include "qv/equiv_quiver_matrix.h"
#include "qv/equivalence_checker.h"
#include "qv/template_exchange_graph.h"

#include <armadillo>

namespace refl {
namespace cartan_exchange {
struct CartanQuiver {
  cluster::EquivQuiverMatrix quiver;
  arma::Mat<int> cartan;
  bool fully_compatible;
};
struct Mutator {
  void operator()(CartanQuiver const* const initial,
                  size_t k,
                  CartanQuiver& output) const;
};
struct Equiv {
  bool operator()(CartanQuiver const* const lhs,
                  CartanQuiver const* const rhs) const;
};
struct Equal {
  bool operator()(CartanQuiver const* const lhs,
                  CartanQuiver const* const rhs) const;
};
struct Hash {
  size_t operator()(CartanQuiver const* const quiver) const;
};
struct NewInstance {
  CartanQuiver* operator()(size_t const size) const;
};
struct DontMutateNonCompatible {
  bool operator()(CartanQuiver const* const mptr, int /*vertex*/) const;
};
using CartanVertex = cluster::exchange_graph::detail::
    VertexInfo<CartanQuiver, Hash, Equiv, Equal, Mutator, NewInstance>;
using CartanGraphInfo = cluster::exchange_graph::detail::GraphInfo<
    DontMutateNonCompatible,
    cluster::exchange_graph::detail::NeverStop>;
}
typedef cluster::exchange_graph::Graph<cartan_exchange::CartanVertex,
                                       cartan_exchange::CartanGraphInfo>
    CartanExchangeGraph;
}

#endif
