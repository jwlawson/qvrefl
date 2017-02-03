#include "cartan_exchange_graph.h"
#include "compatible_cartan_iterator.h"

int main() {
  std::string quiver_str =
		//"{ { 0 1 1 -1 0 } { -1 0 -1 0 1 } { -1 1 0 1 -1 } { 1 0 -1 0 -1 } { 0 -1 1 1 0 } }";
		"{ { 0 -1 -1 1 0 } { 1 0 -1 -1 1 } { 1 1 0 0 -1 } { -1 1 0 0 -1 } { 0 -1 1 1 0 } }";
		/*
      "{ { 0 1 1 0 -1 -1 }"
      "  { -1 0 -1 1 0 1 }"
      "  { -1 1 0 -1 1 0 }"
      "  { 0 -1 1 0 1 -1 }"
      "  { 1 0 -1 -1 0 1 }"
      "  { 1 -1 0 1 -1 0 } }";
			*/
	cluster::EquivQuiverMatrix m(quiver_str);
	refl::CompatibleCartanIterator init_cartan_iter(m);

	if(!init_cartan_iter.has_next()) {
		std::cout << "No compatible cartan found\n";
		return 1;
	}
	refl::cartan_exchange::CartanQuiver initial { m, init_cartan_iter.next(), true };

	refl::CartanExchangeGraph graph(initial, m.num_rows());

	for(auto pair : graph) {
		refl::cartan_exchange::CartanQuiver const& quiver = *pair.first;
		std::cout << quiver.quiver << '\n';
		quiver.cartan.print();
		std::cout << std::boolalpha << quiver.fully_compatible << "\n\n";
	}
}
