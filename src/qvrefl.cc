/**
 * qvrefl.cc
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
#include <unistd.h>
#include <iostream>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include "qv/equiv_mutation_class_loader.h"
#include "qv/mutation_class_loader.h"
#include "qv/stream_iterator.h"

#include "cartan_equiv.h"
#include "cartan_iterator.h"
#include "cartan_mutator.h"
#include "compatible_cartan.h"
#include "filtered_iterator.h"
#include "mutation_star.h"
#include "semi_positive_filter.h"
#include "unique_matrix_filter.h"
#include "util.h"
#include "vector_mutator.h"

namespace refl {
namespace {
using MatrixVec = std::vector<arma::Mat<int>>;
using UniqueCartanIter =
    FilteredIterator<CartanIterator, arma::Mat<int>, UniqueMatrixFilter>;
using SemiPosIter =
    FilteredIterator<UniqueCartanIter, arma::Mat<int>, SemiPositiveFilter>;

auto
get_cartan_iterator(cluster::QuiverMatrix const& q) {
  CartanIterator cartan(q);
  // UniqueCartanIter unique(std::move(cartan));
  // return SemiPosIter(std::move(unique));
  return UniqueCartanIter(std::move(cartan));
}
/**
 * For each vertex in the given quiver, find a semipositive quasi-Cartan which
 * serves that vertex. These matrices are returned in a vector indexed by the
 * label of the vertices.
 *
 * If no such quasi-Cartan companion exists for a given vertex, then there will
 * still be a matrix in the vector, however all its values will be zero.
 */
MatrixVec
find_serving_cartans(cluster::QuiverMatrix const& q) {
  MatrixVec result(q.num_cols(), arma::Mat<int>(q.num_rows(), q.num_cols()));
  boost::dynamic_bitset<> found{q.num_cols(), 0};
  const MutationStar star(q);
  const arma::mat initial_vecs(q.num_rows(), q.num_cols(), arma::fill::eye);
  arma::mat mutated_vecs(q.num_rows(), q.num_cols());
  arma::mat gram_matrix(q.num_rows(), q.num_cols());

  CompatibleCartan compatible;

  auto cartan_iter = SemiPosIter(get_cartan_iterator(q));

  while (cartan_iter.has_next() && !found.all()) {
    arma::Mat<int> const& AQ = cartan_iter.next();
    VectorMutator vmut(q, AQ);
    int_fast16_t ncols = AQ.n_cols;
    for (int_fast16_t i = 0; i < ncols; ++i) {
      vmut.mutate(initial_vecs, i, mutated_vecs);
      gram_matrix = util::gram(mutated_vecs, AQ);
      if (!found[i] && compatible(star.qv(i), gram_matrix)) {
        result[i] = AQ;
        found[i] = true;
      }
    }
  }
  return result;
}
/**
 * Print the first semipositive quasi-Cartan companion serving each vertex of
 * the given quiver. If no such companion exists, then print that the vertex is
 * not served by a Cartan matrix.
 */
void
output_serving_cartans(cluster::QuiverMatrix const& q,
                       std::ostream& os = std::cout) {
  auto res = refl::find_serving_cartans(q);
  int_fast16_t size = res.size();
  for (int_fast16_t i = 0; i < size; ++i) {
    if (res[i].at(0, 0) != 2) {
      os << i << " is not served by any semi-positive quasi-Cartan"
         << os.widen('\n');
    } else {
      std::string label;
      label.append(std::to_string(i)).append(" served by");
      res[i].print(os, label);
    }
  }
}
/**
 * Find the first semipositive quasi-Cartan fully compatible with q.
 *
 * If no such quasi-Cartan exists then the first item in the returned pair is
 * false, otherwise the first item is true and the second item is this first
 * Cartan matrix.
 */
std::pair<bool, arma::Mat<int>>
first_compatible_cartan(cluster::QuiverMatrix const& q) {
  bool result = false;
  arma::Mat<int> cartan;
  const MutationStar star(q);
  arma::Mat<int> gram_matrix(q.num_rows(), q.num_cols());

  CompatibleCartan compatible;

  auto cartan_iter = get_cartan_iterator(q);
  SemiPositiveFilter semipos;

  CartanMutator cmut(q);
  while (cartan_iter.has_next() && !result) {
    bool is_comp = true;
    arma::Mat<int> const& AQ = cartan_iter.next();
    int_fast16_t ncols = AQ.n_cols;
    for (int_fast16_t i = 0; is_comp && i < ncols; ++i) {
      cmut(AQ, i, gram_matrix);
      is_comp = compatible(star.qv(i), gram_matrix);
    }
    result = is_comp && semipos(AQ);
    if (result) {
      cartan = AQ;
    }
  }
  return {result, cartan};
}
/* Check if a semipositive cartan matrix is fully compatible with the quiver. */
bool
check_compatible(cluster::QuiverMatrix const& q, arma::Mat<int> const& AQ,
                 std::ostream& os = std::cout) {
  bool result = true;
  const MutationStar star(q);
  arma::Mat<int> gram_matrix(q.num_rows(), q.num_cols());

  CompatibleCartan compatible;
  if (!compatible(q, AQ)) {
    os << "Initial matrix not compatible" << os.widen('\n');
    result = false;
  }
  CartanMutator cmut(q);
  int_fast16_t ncols = AQ.n_cols;
  for (int_fast16_t i = 0; i < ncols; ++i) {
    cmut(AQ, i, gram_matrix);
    if (!compatible(star.qv(i), gram_matrix)) {
      os << "Mutation at " << i << " not compatible" << os.widen('\n');
      result = false;
    }
  }
  return result;
}
/**
 * Check if the semipositive quasi-Cartan matrix AQ is fully compatible with the
 * given quiver.
 *
 * To check if it is, the Cartan is mutated in every possible direction, and for
 * each mutation is checked to be a companion of the corresponding mutated
 * quiver.
 */
bool
repeat_check_compatible(cluster::QuiverMatrix const& q,
                        MutationStar const& star, arma::Mat<int> const& AQ) {
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
/**
 * Check whether the provided quasi-Cartan matrix is a companion of the given
 * quiver.
 *
 * The Cartan is provided as a QuiverMatrix, so will be converted into an arma
 * matrix. This will involve some allocation.
 */
bool
check_compatible(cluster::QuiverMatrix const& q,
                 cluster::QuiverMatrix const& cartan,
                 std::ostream& os = std::cout) {
  arma::Mat<int> const AQ = util::to_arma(cartan);
  bool result = check_compatible(q, AQ, os);
  return result;
}
/**
 * Compute all fully-compatible, semipositive quasi-Cartan companions of the
 * given quiver. These Cartan matrices are returned in a vector.
 *
 * Only distinct, i.e. non-equivalent, Cartan matrices will be returned. (Where
 * equivalence is given by flipping signs at a vertex.)
 */
std::vector<arma::Mat<int>>
all_compatible(cluster::QuiverMatrix const& q) {
  auto cartan_iter = get_cartan_iterator(q);
  MutationStar const star(q);
  CartanEquiv equiv;
  SemiPositiveFilter semipos;

  std::vector<arma::Mat<int>> result;
  result.reserve(2);

  do {
    auto const& first = cartan_iter.next();
    if (repeat_check_compatible(q, star, first) && semipos(first)) {
      result.push_back(first);
    }
  } while (result.empty() && cartan_iter.has_next());

  if (!result.empty()) {
    while (cartan_iter.has_next()) {
      auto const& n = cartan_iter.next();
      if (repeat_check_compatible(q, star, n) && semipos(n) &&
          std::find_if(result.begin(), result.end(),
                       [&n, &equiv](arma::Mat<int> const& c) {
                         return equiv(c, n);
                       }) == result.end()) {
        result.push_back(n);
      }
    }
  }
  return result;
}
/**
 * Check whether all fully-compatible, semi-positive quasi-Cartan companions of
 * the given quiver are equivalent (up to flipping signs at vertices). If not
 * all the possble companions are printed to the supplied ostream.
 */
bool
check_all_compatible_equiv(cluster::QuiverMatrix const& q,
                           std::ostream& os = std::cout) {
  auto compatible_cartans = all_compatible(q);
  bool all_equiv = compatible_cartans.size() == 1;

  if (!all_equiv) {
    for (auto const& n : compatible_cartans) {
      n.print(os, "Compatible:");
    }
  }
  return all_equiv;
}
}  // namespace
}  // namespace refl
enum Func { ListCompatible, SameCompatible, CheckSingle, CompatibleEquiv };
void
usage() {
  std::cout << "qvrefl -cels [-m matrix] [-i in_file] [-a cartan]"
            << std::cout.widen('\n');
  std::cout << "  -m Specify matrix to find quasi-Cartan companions of"
            << std::cout.widen('\n');
  std::cout << "  -i Specify input file of matrices to read"
            << std::cout.widen('\n');
  std::cout << "  -c Check the whole mutation class of the matrix (only with "
               "-l or -s)"
            << std::cout.widen('\n');
  std::cout << "  -e Check that all fully compatible cartans are equivalent"
            << std::cout.widen('\n');
  std::cout << "  -l List the first quasi-Cartan serving each vertex"
            << std::cout.widen('\n');
  std::cout << "  -s Check if there is a quasi-Cartan serving all vertices"
            << std::cout.widen('\n');
  std::cout
      << "  -a Specify a cartan matrix to check whether it serves every vertex"
      << std::cout.widen('\n');
  std::cout.flush();
}
int
main(int argc, char* argv[]) {
  std::string matrix;
  std::string cartan;
  Func function = Func::ListCompatible;
  std::string input;
  bool mut_class = false;
  int c;
  while ((c = getopt(argc, argv, "elsm:i:hca:")) != -1) {
    switch (c) {
      case 'e':
        function = Func::CompatibleEquiv;
        break;
      case 'm':
        matrix = optarg;
        break;
      case 'a':
        cartan = optarg;
        function = Func::CheckSingle;
        break;
      case 'l':
        function = Func::ListCompatible;
        break;
      case 's':
        function = Func::SameCompatible;
        break;
      case 'i':
        input = optarg;
        break;
      case 'c':
        mut_class = true;
        break;
      case 'h':
      case '?':
      default:
        usage();
        return 1;
    }
  }
  if (matrix.length() < 4 && input.empty()) {
    usage();
    return 2;
  }
  if (function == Func::ListCompatible) {
    if (matrix.length() > 0) {
      if (mut_class) {
        cluster::EquivQuiverMatrix q(matrix);
        cluster::EquivMutationClassLoader cl(q);
        while (cl.has_next()) {
          auto quiver = cl.next_ptr();
          std::cout << *quiver << ":" << std::cout.widen('\n');
          refl::output_serving_cartans(*quiver);
        }
      } else {
        cluster::QuiverMatrix q(matrix);
        refl::output_serving_cartans(q);
      }
    } else {
      std::ifstream inf;
      inf.open(input);
      if (!inf.is_open()) {
        std::cerr << "Could not open file " << input << std::endl;
        return 1;
      }
      cluster::StreamIterator<cluster::QuiverMatrix> iter(inf);
      while (iter.has_next()) {
        auto quiver = iter.next();
        std::cout << *quiver << std::cout.widen('\n');
        refl::output_serving_cartans(*quiver);
      }
    }
  } else if (function == Func::SameCompatible) {
    if (matrix.length() > 0) {
      if (mut_class) {
        cluster::EquivQuiverMatrix q(matrix);
        cluster::EquivMutationClassLoader cl(q);
        while (cl.has_next()) {
          auto quiver = cl.next_ptr();
          bool result = refl::first_compatible_cartan(*quiver).first;
          std::cout << (result ? "True: " : "False: ") << *quiver
                    << std::cout.widen('\n');
        }
      } else {
        cluster::QuiverMatrix q(matrix);
        auto result = refl::first_compatible_cartan(q);
        std::cout << (result.first ? result.second : "False")
                  << std::cout.widen('\n');
      }
    } else {
      std::ifstream inf;
      inf.open(input);
      if (!inf.is_open()) {
        std::cerr << "Could not open file " << input << std::endl;
        return 1;
      }
      cluster::StreamIterator<cluster::QuiverMatrix> iter(inf);
      while (iter.has_next()) {
        auto quiver = iter.next();
        bool result = refl::first_compatible_cartan(*quiver).first;
        std::cout << (result ? "True: " : "False: ") << *quiver
                  << std::cout.widen('\n');
      }
    }
  } else if (function == Func::CheckSingle) {
    if (matrix.empty() || cartan.empty()) {
      usage();
      return 4;
    }
    refl::check_compatible(cluster::QuiverMatrix{matrix},
                           cluster::QuiverMatrix{cartan});
  } else if (function == Func::CompatibleEquiv) {
    if (!matrix.empty()) {
      std::cout << std::boolalpha
                << refl::check_all_compatible_equiv(
                       cluster::QuiverMatrix{matrix})
                << std::cout.widen('\n');
    } else if (!input.empty()) {
      std::ifstream inf;
      inf.open(input);
      if (!inf.is_open()) {
        std::cerr << "Could not open file " << input << std::endl;
        return 1;
      }
      cluster::StreamIterator<cluster::QuiverMatrix> iter(inf);
      while (iter.has_next()) {
        auto quiver = iter.next();
        std::cout << std::boolalpha << refl::check_all_compatible_equiv(*quiver)
                  << ": " << *quiver << std::cout.widen('\n');
      }
    } else {
      usage();
      return 4;
    }
  } else {
    usage();
  }
  std::cout.flush();
  return 0;
}
