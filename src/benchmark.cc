#include <benchmark/benchmark.h>

#include "cartan_equiv.h"
#include "cartan_mutator.h"
#include "vector_mutator.h"

static void
CartanEquivSameCheck(benchmark::State& state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    ::benchmark::DoNotOptimize(equiv(a, a));
  }
}
BENCHMARK(CartanEquivSameCheck);

static void
CartanEquivOneSwitch(benchmark::State& state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> b{
      {2, -1, -1, -1}, {-1, 2, 1, 1}, {-1, 1, 2, 1}, {-1, 1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    ::benchmark::DoNotOptimize(equiv(a, a));
  }
}
BENCHMARK(CartanEquivOneSwitch);

static void
CartanEquivTwoSwitch(benchmark::State& state) {
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  arma::Mat<int> d{
      {2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};
  refl::CartanEquiv equiv;

  while (state.KeepRunning()) {
    ::benchmark::DoNotOptimize(equiv(a, a));
  }
}
BENCHMARK(CartanEquivTwoSwitch);

static void
CartanMutatorA3(benchmark::State& state) {
  cluster::QuiverMatrix a3{"{ { 0 1 -1 } { -1 0 1 } { 1 -1 0 } }"};
  arma::Mat<int> a{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
  refl::CartanMutator mut(a3);

  arma::Mat<int> output;

  while (state.KeepRunning()) {
    mut(a, 0, output);
  }
}
BENCHMARK(CartanMutatorA3);

static void
CartanMutatorA4(benchmark::State& state) {
  cluster::QuiverMatrix a4{
      "{ { 0 1 0 0 } { -1 0 1 0 } { 0 -1 0 1 0 } { 0 0 -1 0 } }"};
  arma::Mat<int> a{{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}};
  refl::CartanMutator mut(a4);

  arma::Mat<int> output;

  while (state.KeepRunning()) {
    mut(a, 0, output);
  }
}
BENCHMARK(CartanMutatorA4);

static void
CartanMutatorA5(benchmark::State& state) {
  cluster::QuiverMatrix a5{
      "{ { 0 1 0 0 0 } { -1 0 1 0 0 } { 0 -1 0 1 0 0 } { "
      "0 0 -1 0 1 } { 0 0 0 -1 0 } }"};
  arma::Mat<int> a{{2, 1, 0, 0, 0},
                   {1, 2, 1, 0, 0},
                   {0, 1, 2, 1, 0},
                   {0, 0, 1, 2, 1},
                   {0, 0, 0, 1, 2}};
  refl::CartanMutator mut(a5);

  arma::Mat<int> output;

  while (state.KeepRunning()) {
    mut(a, 0, output);
  }
}
BENCHMARK(CartanMutatorA5);

static void
VectorMutateA3(benchmark::State& state) {
  arma::Mat<int> a = {{2, 1, 0}, {1, 2, 1}, {0, 1, 2}};
  cluster::QuiverMatrix q("{ { 0 1 0 } { -1 0 1 } { 0 -1 0 } }");
  refl::VectorMutator vm(q, a);
  arma::Mat<int> vecs = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  arma::Mat<int> res(3, 3);
  while (state.KeepRunning()) {
    vm.mutate(vecs, 0, res);
  }
}
BENCHMARK(VectorMutateA3);

static void
VectorMutateA4(benchmark::State& state) {
  arma::Mat<int> a = {{2, 1, 0, 0}, {1, 2, 1, 0}, {0, 1, 2, 1}, {0, 0, 1, 2}};
  cluster::QuiverMatrix q(
      "{ { 0 1 0 0 } { -1 0 1 0 } { 0 -1 0 1 } { 0 0 -1 0 } }");
  refl::VectorMutator vm(q, a);
  arma::Mat<int> vecs = {
      {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  arma::Mat<int> res(4, 4);
  while (state.KeepRunning()) {
    vm.mutate(vecs, 0, res);
  }
}
BENCHMARK(VectorMutateA4);

static void
VectorMutateA5(benchmark::State& state) {
  arma::Mat<int> a = {{2, 1, 0, 0, 0},
                      {1, 2, 1, 0, 0},
                      {0, 1, 2, 1, 0},
                      {0, 0, 1, 2, 1},
                      {0, 0, 0, 1, 2}};
  cluster::QuiverMatrix q(
      "{ { 0 1 0 0 0 } { -1 0 1 0 0 } { 0 -1 0 1 0 } { 0 0 "
      "-1 0 1 } { 0 0 0 -1 0 } }");
  refl::VectorMutator vm(q, a);
  arma::Mat<int> vecs = {{1, 0, 0, 0, 0},
                         {0, 1, 0, 0, 0},
                         {0, 0, 1, 0, 0},
                         {0, 0, 0, 1, 0},
                         {0, 0, 0, 0, 1}};
  arma::Mat<int> res(5, 5);
  while (state.KeepRunning()) {
    vm.mutate(vecs, 0, res);
  }
}
BENCHMARK(VectorMutateA5);

static void
VectorMutateA5Double(benchmark::State& state) {
  arma::Mat<int> a = {{2, 1, 0, 0, 0},
                      {1, 2, 1, 0, 0},
                      {0, 1, 2, 1, 0},
                      {0, 0, 1, 2, 1},
                      {0, 0, 0, 1, 2}};
  cluster::QuiverMatrix q(
      "{ { 0 1 0 0 0 } { -1 0 1 0 0 } { 0 -1 0 1 0 } { 0 0 "
      "-1 0 1 } { 0 0 0 -1 0 } }");
  refl::VectorMutator vm(q, a);
  arma::Mat<double> vecs = {{1, 0, 0, 0, 0},
                            {0, 1, 0, 0, 0},
                            {0, 0, 1, 0, 0},
                            {0, 0, 0, 1, 0},
                            {0, 0, 0, 0, 1}};
  arma::Mat<double> res(5, 5);
  while (state.KeepRunning()) {
    vm.mutate(vecs, 0, res);
  }
}
BENCHMARK(VectorMutateA5Double);
BENCHMARK_MAIN();
