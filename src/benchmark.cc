#include <benchmark/benchmark.h>

#include "cartan_equiv.h"

static void CartanEquivSameCheck(benchmark::State &state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    volatile bool test = equiv(a, a);
  }
}
BENCHMARK(CartanEquivSameCheck);

static void CartanEquivOneSwitch(benchmark::State &state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{
      {2, -1, -1, -1}, {-1, 2, 1, 1}, {-1, 1, 2, 1}, {-1, 1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    volatile bool test = equiv(a, b);
  }
}
BENCHMARK(CartanEquivOneSwitch);

static void CartanEquivTwoSwitch(benchmark::State &state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> d{
      {2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    volatile bool test = equiv(a, d);
  }
}
BENCHMARK(CartanEquivTwoSwitch);

BENCHMARK_MAIN();
