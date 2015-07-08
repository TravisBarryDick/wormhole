#ifndef RPTS_H
#define RPTS_H

#include <memory>
#include <random>
#include <vector>

#include "kmeans_helper.h"
#include "local/row_block.h"
#include "local/data.h"

namespace dddml {

namespace rpt_impl { // Namespace for the implementation

using namespace std;
using namespace dmlc;
using namespace dmlc::data;

enum RouteDirection { LEFT, RIGHT };

////////////////////
// RPTSplit class //
////////////////////

class RPTSplit {
public:
  /* Constructs a placeholder RPTSplit used in leaf-nodes. */
  RPTSplit();

  /* Randomly generates a new RPTSplit. The split point t is chosen based on the
   * rows of data given by idxs. */
  template <typename IndexType>
  RPTSplit(mt19937_64 &rng, size_t dimension, const RowBlock<IndexType> data,
           const vector<IndexType> &idxs);

  /* Determines which direction a given row should be routed by this
   * split */
  template <typename IndexType> RouteDirection route(Row<IndexType> row);

private:
  vector<real_t> d; // The direction for the split
  real_t t;         // The split point
};

// *** RPTSplit Definitions *** //

RPTSplit::RPTSplit() : d(0), t(0) {}

template <typename IndexType>
RPTSplit::RPTSplit(mt19937_64 &rng, size_t dimension,
                   const RowBlock<IndexType> data,
                   const vector<IndexType> &idxs)
    : d(dimension) {
  // sample a random direction d
  auto std_normal = normal_distribution<real_t>();
  for (size_t i = 0; i < dimension; i++) {
    d[i] = std_normal(rng);
  }
  // project each sample onto d
  vector<real_t> ps(idxs.size());
  for (int i = 0; i < idxs.size(); ++i) {
    ps[i] = data[idxs[i]].SDot(d.data(), dimension);
  }
  // Pick a random fractile between 1/4 and 3/4
  real_t fractile = std::uniform_real_distribution<real_t>(0.25, 0.75)(rng);
  // Use std::nth_element to find the split threshold
  int findex = int(ps.size() * fractile);
  nth_element(ps.begin(), ps.begin() + findex, ps.end());
  t = ps[findex];
};

template <typename IndexType>
inline auto RPTSplit::route(Row<IndexType> row) -> RouteDirection {
  real_t p = row.SDot(d.data(), d.size());
  return (p > t) ? RIGHT : LEFT;
}

///////////////////
// RPTNode class //
///////////////////

template <typename IndexType> class RPTNode {
public:
  /* Constructs a leaf node containing the given indices. */
  RPTNode(const vector<IndexType> &idxs);

  /* Constructs an internal node with the given split and given subtrees. */
  RPTNode(const RPTSplit &split, unique_ptr<RPTNode> &left,
          unique_ptr<RPTNode> &right);

  /* Returns the depth of the tree */
  auto depth() -> size_t;

  /* Returns the row indices in the same leaf as row. */
  auto search(Row<IndexType> row) -> const vector<IndexType> *;

private:
  unique_ptr<RPTNode> left;
  unique_ptr<RPTNode> right;
  RPTSplit split;
  vector<IndexType> idxs;
};

// *** RPTNode Definitions *** //

template <typename IndexType>
RPTNode<IndexType>::RPTNode(const vector<IndexType> &idxs)
    : idxs(idxs) {}

template <typename IndexType>
RPTNode<IndexType>::RPTNode(const RPTSplit &split, unique_ptr<RPTNode> &left,
                            unique_ptr<RPTNode> &right)
    : left(std::move(left)), right(std::move(right)), split(split) {}

template <typename IndexType> inline size_t RPTNode<IndexType>::depth() {
  return !(left && right) ? 1 : max(left->depth(), right->depth());
}

template <typename IndexType>
inline const vector<IndexType> *RPTNode<IndexType>::search(Row<IndexType> row) {
  if (!(left && right))
    return &idxs;
  else
    return (split.route(row) == LEFT ? left : right)->search(row);
};

//////////////////////////////////////////
// Helper functions for making RPTNodes //
//////////////////////////////////////////

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data,
                                           vector<IndexType> &idxs) {

  if (idxs.size() <= n0) {
    return unique_ptr<RPTNode<IndexType>>(new RPTNode<IndexType>(idxs));
  } else {
    auto split = RPTSplit(rng, dimension, data, idxs);
    vector<IndexType> left_idxs(0);
    vector<IndexType> right_idxs(0);
    for (IndexType idx : idxs) {
      if (split.route(data[idx]) == LEFT)
        left_idxs.push_back(idx);
      else
        right_idxs.push_back(idx);
    }
    auto left = make_rptree(rng, dimension, n0, data, left_idxs);
    auto right = make_rptree(rng, dimension, n0, data, right_idxs);
    return unique_ptr<RPTNode<IndexType>>(
        new RPTNode<IndexType>(split, left, right));
  }
}

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data) {
  vector<IndexType> idxs(data.size);
  for (size_t i = 0; i < data.size; ++i)
    idxs[i] = static_cast<IndexType>(i);
  return make_rptree(rng, dimension, n0, data, idxs);
}
};

///////////////////////////////
// RandomPartitionTree class //
///////////////////////////////

template <typename IndexType> class RandomPartitionTree {
public:
  RandomPartitionTree(std::mt19937_64 &rng, int dimension, int n0,
                      dmlc::data::RowBlockContainer<IndexType> &data);

  RandomPartitionTree(std::mt19937_64 &rng, int dimension, int n0,
                      dmlc::data::RowBlockContainer<IndexType> &data,
                      std::vector<IndexType> &idxs);

  auto depth() -> size_t;

  auto find_nn(dmlc::Row<IndexType> row) -> IndexType;

  auto get_rowblock() -> const dmlc::data::RowBlockContainer<IndexType> &;

private:
  dmlc::data::RowBlockContainer<IndexType> data;
  std::unique_ptr<rpt_impl::RPTNode<IndexType>> tree;
};

// *** RandomPartitionTree Definitions *** //

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(
    std::mt19937_64 &rng, int dimension, int n0,
    dmlc::data::RowBlockContainer<IndexType> &data)
    : data(data),
      tree(rpt_impl::make_rptree(rng, dimension, n0, data.GetBlock())) {}

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(
    std::mt19937_64 &rng, int dimension, int n0,
    dmlc::data::RowBlockContainer<IndexType> &data,
    std::vector<IndexType> &idxs)
    : data(data),
      tree(rpt_impl::make_rptree(rng, dimension, n0, data.GetBlock(), idxs)) {}

template <typename IndexType> size_t RandomPartitionTree<IndexType>::depth() {
  return tree->depth();
}

template <typename IndexType>
IndexType RandomPartitionTree<IndexType>::find_nn(dmlc::Row<IndexType> row) {
  auto leaf_idxs = tree->search(row);
  real_t min_dist = std::numeric_limits<real_t>::infinity();
  IndexType best_idx;
  for (IndexType idx : *leaf_idxs) {
    real_t dist = squareDist(row, data.GetBlock()[idx]);
    if (dist < min_dist) {
      min_dist = dist;
      best_idx = idx;
    }
  }
  return best_idx;
}

template <typename IndexType>
const dmlc::data::RowBlockContainer<IndexType> &
RandomPartitionTree<IndexType>::get_rowblock() {
  return data;
}
};
#endif
