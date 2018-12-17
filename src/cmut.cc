#include "cartan_mutator.h"
#include "compatible_cartan.h"
#include "compatible_cartan_iterator.h"
#include "permuted_cartan_equiv.h"

#include "qv/equiv_quiver_matrix.h"

#include <iostream>
#include <unordered_map>
#include <vector>

#include "armadillo"

struct CartanInfo {
  int index;
  std::vector<arma::Mat<int>> cartans;
};

typedef cluster::EquivQuiverMatrix Quiver;
typedef std::unordered_map<Quiver, CartanInfo> Map;

struct State {
  int index;
  Quiver const* quiver;
  arma::Mat<int> const* cartan;
};
State
insert_quiver(Map& map, Quiver q, CartanInfo info) {
  auto ins = map.emplace(std::move(q), std::move(info));
  auto& first_iter = ins.first;
  auto& c_info = ins.first->second;
  auto& cartan_vec = c_info.cartans;
  size_t num_cartans = cartan_vec.size();

  return {c_info.index, &(first_iter->first), &(cartan_vec[num_cartans - 1])};
}
State
insert_quiver(Map& map, int index, Quiver q, arma::Mat<int> cartan) {
  CartanInfo info{index, {std::move(cartan)}};
  return insert_quiver(map, q, std::move(info));
}

int
main(int argc, char* argv[]) {
  Map map;
  std::string quiver_str;
  if (argc > 1) {
    quiver_str = argv[1];
  } else {
    quiver_str =
        "{ { 0 1 1 0 -1 -1 }"
        "  { -1 0 -1 1 0 1 }"
        "  { -1 1 0 -1 1 0 }"
        "  { 0 -1 1 0 1 -1 }"
        "  { 1 0 -1 -1 0 1 }"
        "  { 1 -1 0 1 -1 0 } }";
  }
  Quiver initial(quiver_str);
  refl::CompatibleCartanIterator cartans(initial);
  CartanInfo initial_info{0, {}};

  if (!cartans.has_next()) {
    std::cout << "Quiver has no fully compatible cartans!\n";
    return 1;
  }
  while (cartans.has_next()) {
    initial_info.cartans.push_back(cartans.next());
  }
  State current =
      insert_quiver(map, std::move(initial), std::move(initial_info));

  int index = 0;
  bool should_continue = true;
  std::string inp;
  while (should_continue) {
    std::cin >> inp;
    bool mutate = false;
    bool print = false;
    int mut;
    switch (inp[0]) {
      case 'q':
        should_continue = false;
        break;
      case 'r':
        break;
      case 'n':
        break;
      case 's':
        print = true;
        break;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        mutate = true;
        mut = std::stoi(inp);
        break;
    }
    if (mutate && mut < current.quiver->num_rows()) {
      std::cout << "Mutating in " << mut << '\n';
      Quiver new_quiver(current.quiver->num_rows(), current.quiver->num_cols());
      arma::Mat<int> new_cartan(current.quiver->num_rows(),
                                current.quiver->num_cols());
      current.quiver->mutate(mut, new_quiver);
      refl::CartanMutator cmut(*current.quiver);
      cmut(*current.cartan, mut, new_cartan);

      State next;
      {
        auto pos = map.find(new_quiver);
        if (pos != map.end()) {
          std::cout << "Mutated quiver same as: " << pos->second.index << '\n';
          next.index = pos->second.index;
          next.quiver = std::addressof(pos->first);
          refl::PermutedCartanEquiv equiv;
          auto& cartan_vec = pos->second.cartans;
          auto permutation = new_quiver.get_permutation(pos->first);
          auto cequiv = [&new_cartan, &equiv,
                         &permutation](arma::Mat<int> const& other) {
            return equiv(new_cartan, other, permutation);
          };
          {
            auto cartan_pos =
                std::find_if(cartan_vec.begin(), cartan_vec.end(), cequiv);
            if (cartan_pos != cartan_vec.end()) {
              std::cout << "Cartan equivalent to: "
                        << std::distance(cartan_vec.begin(), cartan_pos)
                        << '\n';
              next.cartan = std::addressof(*cartan_pos);
            } else {
              std::cout << "New cartan matrix\n";
              refl::FullyCompatibleCheck compatible;
              if (!compatible(new_quiver, new_cartan)) {
                std::cout << "Cartan not fully compatible!\n";
              }
              cartan_vec.push_back(std::move(new_cartan));
              next.cartan = std::addressof(cartan_vec[cartan_vec.size() - 1]);
            }
          }
        } else {
          ++index;
          std::cout << "New quiver. Index: " << index << '\n';
          refl::FullyCompatibleCheck compatible;
          if (!compatible(new_quiver, new_cartan)) {
            std::cout << "Cartan not fully compatible!\n";
          }
          next = insert_quiver(map, index, std::move(new_quiver),
                               std::move(new_cartan));
        }
      }
      current = std::move(next);
    }
    if (print) {
      std::cout << "Index: " << current.index << '\n';
      std::cout << *current.quiver << '\n';
      current.cartan->print("Cartan:");
    }
  }
}
