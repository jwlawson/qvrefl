#include "cartan_exchange_graph.h"

#include "compatible_cartan_iterator.h"

#include "gtest/gtest.h"

namespace refl {
TEST(CartanExchangeGraph, NotFullyCompatible) {
  std::string str =
      "{ { 0 -1 -1 1 0 } { 1 0 -1 -1 1 } { 1 1 0 0 -1 } { -1 1 0 0 -1 } { 0 -1 "
      "1 1 0 } }";
  cluster::EquivQuiverMatrix m(str);
  refl::CompatibleCartanIterator init_cartan_iter(m);

  ASSERT_TRUE(init_cartan_iter.has_next());

  refl::cartan_exchange::CartanQuiver initial{m, init_cartan_iter.next(), true};

  refl::CartanExchangeGraph graph(initial, m.num_rows());

  int n_verts = 0;
  int n_comp = 0;
  for (auto pair : graph) {
    n_verts++;
    refl::cartan_exchange::CartanQuiver const& quiver = *pair.first;
    if (quiver.fully_compatible) {
      n_comp++;
    }
  }
  EXPECT_EQ(12, n_verts);
  EXPECT_EQ(8, n_comp);
}
TEST(CartanExchangeGraph, EquivCheck) {
  std::string str =
      "{ { 0 1 0 0 0 } { -1 0 0 0 -1 } {  0 0 0 1 0 } { 0 0 -1 0 1 } { 0 1 0 "
      "-1 0 } }";

  cluster::EquivQuiverMatrix m(str);
  arma::Mat<int> cartan = {2, 1, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 2,
                           1, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 2};
  cartan.reshape(5, 5);

  refl::cartan_exchange::CartanQuiver quiv1{m, cartan, true};
  refl::cartan_exchange::CartanQuiver quiv2{m, cartan, true};

  refl::cartan_exchange::Hash hash;
  size_t hash1 = hash(&quiv1);
  size_t hash2 = hash(&quiv2);
  EXPECT_EQ(hash1, hash2);

  refl::cartan_exchange::Equiv equiv;
  EXPECT_TRUE(equiv(&quiv1, &quiv2));
}
}  // namespace refl
