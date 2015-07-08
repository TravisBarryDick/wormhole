#ifndef RPTS_H
#define RPTS_H

#include "base.h"

#include <memory>
#include <random>
#include <vector>
#include "kmeans_helper.h"

#if DISTRIBUTED

#include "data/row_block.h"
#include "dmlc/data.h"
#include "dmlc/io.h"

#else

#include "local/row_block.h"
#include "local/data.h"

#endif
namespace dddml {

namespace rpt_impl { // Namespace for the implementation

using namespace std;
using namespace dmlc;
using namespace dmlc::data;

enum RouteDirection { LEFT, RIGHT };

////////////////////
// RPTSplit class //
////////////////////

template <typename IndexType>
class RPTSplit {
public:
  /* Constructs a placeholder RPTSplit used in leaf-nodes. */
  RPTSplit(std::vector<IndexType> * feature_dict_ptr = NULL) ;

  /* Randomly generates a new RPTSplit. The split point t is chosen based on the
   * rows of data given by idxs. */
  RPTSplit(mt19937_64 &rng, size_t dimension, const RowBlock<IndexType> data,
           const vector<size_t> &idxs, vector<IndexType> *feature_dict_ptr);

  /* Determines which direction a given row should be routed by this
   * split */
  RouteDirection route(Row<IndexType> row);
  
  /* Disk IO functions */
  void Save(dmlc::Stream *fo);
  static unique_ptr<RPTSplit> Load(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr);
  
private:
  vector<real_t> d; // The direction for the split
  real_t t;         // The split point
  vector<IndexType> *feature_dict_ptr; // feature map
};

// *** RPTSplit Definitions *** //

template <typename IndexType>
RPTSplit<IndexType>::RPTSplit(std::vector<IndexType> *feature_dict_ptr) : d(0), t(0), feature_dict_ptr(feature_dict_ptr) {}

template <typename IndexType>
RPTSplit<IndexType>::RPTSplit(mt19937_64 &rng, size_t dimension,
                   const RowBlock<IndexType> data,
                   const vector<size_t> &idxs,
                   vector<IndexType> *feature_dict_ptr)
    : d(dimension) , feature_dict_ptr(feature_dict_ptr) {
  // sample a random direction d
  auto std_normal = normal_distribution<real_t>();
  for (size_t i = 0; i < dimension; i++) {
    d[i] = std_normal(rng);
  }
  // project each sample onto d
  vector<real_t> ps(idxs.size());
  for (size_t i = 0; i < idxs.size(); ++i) {
    ps[i] = data[idxs[i]].SDot(d.data(), dimension); //TODO: Change
  }
  // Pick a random fractile between 1/4 and 3/4
  real_t fractile = std::uniform_real_distribution<real_t>(0.25, 0.75)(rng);
  // Use std::nth_element to find the split threshold
  int findex = int(ps.size() * fractile);
  nth_element(ps.begin(), ps.begin() + findex, ps.end());
  t = ps[findex];
};

template <typename IndexType>
inline auto RPTSplit<IndexType>::route(Row<IndexType> row) -> RouteDirection {
  real_t p = row.SDot(d.data(), d.size()); //TODO: change
  return (p > t) ? RIGHT : LEFT;
}


/* IO functions */
template <typename IndexType>
void RPTSplit<IndexType>::Save(dmlc::Stream *fo){
  //save split point first, and then the split direction
  fo->Write(&t, sizeof(real_t));
  fo->Write(d);
}

template <typename IndexType>
inline unique_ptr<RPTSplit<IndexType>> RPTSplit<IndexType>::Load(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr1){
  //read split point, and then split direction
  auto ptr = unique_ptr<RPTSplit<IndexType>>(new RPTSplit(feature_dict_ptr1));
  fi->Read(&(ptr->t), sizeof(real_t));
  fi->Read(&(ptr->d));
  return ptr;
}


///////////////////
// RPTNode class //
///////////////////

template <typename IndexType> class RPTNode {
public:
  /* Constructs a leaf node containing the given indices. */
  RPTNode(const vector<size_t> &idxs);

  /* Constructs an internal node with the given split and given subtrees. */
  RPTNode(const RPTSplit<IndexType> &split, unique_ptr<RPTNode<IndexType>> &left,
          unique_ptr<RPTNode<IndexType>> &right);

  /* Returns the depth of the tree */
  auto depth() -> size_t;

  /* Returns the row indices in the same leaf as row. */
  auto search(Row<IndexType> row) -> const vector<size_t> *;
  
  /* Disk IO functions */
  void Save(dmlc::Stream *fo);
  static unique_ptr<RPTNode<IndexType>> Load(dmlc::Stream *fi);
  
private:
  unique_ptr<RPTNode<IndexType>> left;
  unique_ptr<RPTNode<IndexType>> right;
  RPTSplit<IndexType> split;
  vector<size_t> idxs;
};

// *** RPTNode Definitions *** //

template <typename IndexType>
RPTNode<IndexType>::RPTNode(const vector<size_t> &idxs)
    : idxs(idxs) {}

template <typename IndexType>
RPTNode<IndexType>::RPTNode(const RPTSplit<IndexType> &split, unique_ptr<RPTNode> &left,
                            unique_ptr<RPTNode> &right)
    : left(std::move(left)), right(std::move(right)), split(split) {}

template <typename IndexType> inline size_t RPTNode<IndexType>::depth() {
  return !(left && right) ? 1 : (1 + max(left->depth(), right->depth()));
}

template <typename IndexType>
inline const vector<size_t> *RPTNode<IndexType>::search(Row<IndexType> row) {
  if (!(left && right))
    return &idxs;
  else
    return (split.route(row) == LEFT ? left : right)->search(row);
};

template <typename IndexType>
inline void RPTNode<IndexType>::Save(dmlc::Stream *fo)
{
  //save ifLeaf = 0 for internal nodes, and ifLeaf = 1 for leaves
  bool isLeaf = (idxs.size() > 0);
  fo->Write(&isLeaf, sizeof(bool));
  
  if (!isLeaf) {
    //save children recursively
    left->Save(fo);
    right->Save(fo);
    //only save split for internal nodes
    split.Save(fo); 
  }
  else fo->Write(idxs); //only save idxs for leaves
}

template <typename IndexType>
static unique_ptr<RPTNode<IndexType>> Load(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr)
{
  bool isLeaf;
  fi->Read(&(isLeaf), sizeof(bool));
  
  if (!isLeaf){
    // read left and right children recursively
    unique_ptr<RPTNode<IndexType>> left = Load(fi, feature_dict_ptr);
    unique_ptr<RPTNode<IndexType>> right = Load(fi, feature_dict_ptr);
    // read split
    auto split_ptr = RPTSplit<IndexType>::Load(fi, feature_dict_ptr);
    return unique_ptr<RPTNode<IndexType>>(new RPTNode<IndexType>(*split_ptr, left, right));
  }
  else{
    //read idxs
    vector<size_t> idxs;
    fi->Read(&idxs);
    return unique_ptr<RPTNode<IndexType>>(new RPTNode<IndexType>(idxs));
  }
}


//////////////////////////////////////////
// Helper functions for making RPTNodes //
//////////////////////////////////////////

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data,
                                           vector<size_t> &idxs,
                                           vector<IndexType> *feature_dict_ptr) {

  if (idxs.size() <= n0) {
    return unique_ptr<RPTNode<IndexType>>(new RPTNode<IndexType>(idxs));
  } else {
    auto split = RPTSplit<IndexType>(rng, dimension, data, idxs, feature_dict_ptr);
    vector<size_t> left_idxs(0);
    vector<size_t> right_idxs(0);
    for (IndexType idx : idxs) {
      if (split.route(data[idx]) == LEFT)
        left_idxs.push_back(idx);
      else
        right_idxs.push_back(idx);
    }
    auto left = make_rptree(rng, dimension, n0, data, left_idxs, feature_dict_ptr);
    auto right = make_rptree(rng, dimension, n0, data, right_idxs, feature_dict_ptr);
    return unique_ptr<RPTNode<IndexType>>(
        new RPTNode<IndexType>(split, left, right));
  }
}

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data, vector<IndexType> *feature_dict_ptr) {
  vector<size_t> idxs(data.size);
  for (size_t i = 0; i < data.size; ++i)
    idxs[i] = i;
  return make_rptree(rng, dimension, n0, data, idxs, feature_dict_ptr);
}


template <typename IndexType>
unique_ptr<RPTNode<IndexType>> read_rptree(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr)
{
  return RPTNode<IndexType>::Load(fi, feature_dict_ptr);
}

}; //namespace rpt_impl

///////////////////////////////
// RandomPartitionTree class //
///////////////////////////////

template <typename IndexType> class RandomPartitionTree {
public:
  RandomPartitionTree(std::mt19937_64 &rng, int dimension, int n0,
                      dmlc::data::RowBlockContainer<IndexType> &data,
                      std::vector<IndexType> &feature_dict);

  RandomPartitionTree(std::mt19937_64 &rng, int dimension, int n0,
                      dmlc::data::RowBlockContainer<IndexType> &data,
                      std::vector<size_t> &idxs,
                      std::vector<IndexType> &feature_dict);

  auto depth() -> size_t;

  auto find_nn(dmlc::Row<IndexType> row) -> IndexType;

  auto get_rowblock() -> const dmlc::data::RowBlockContainer<IndexType> &;
  
  void Save(dmlc::Stream *);
  RandomPartitionTree(dmlc::Stream*, dmlc::Stream*);

private:
  std::vector<IndexType> feature_dict;
  dmlc::data::RowBlockContainer<IndexType> data;
  std::unique_ptr<rpt_impl::RPTNode<IndexType>> tree;
};

// *** RandomPartitionTree Definitions *** //

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(
    std::mt19937_64 &rng, int dimension, int n0,
    dmlc::data::RowBlockContainer<IndexType> &data,
    std::vector<IndexType> &feature_dict)
    : data(data), feature_dict(feature_dict)
{
  tree = std::move(rpt_impl::make_rptree(rng, dimension, n0, data.GetBlock(), &(this->feature_dict)));
}
      

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(
    std::mt19937_64 &rng, int dimension, int n0,
    dmlc::data::RowBlockContainer<IndexType> &data,
    std::vector<size_t> &idxs,
    std::vector<IndexType> &feature_dict)
    : data(data), feature_dict(feature_dict)
{
  tree = rpt_impl::make_rptree(rng, dimension, n0, data.GetBlock(), idxs, &(this->feature_dict));      
}

template <typename IndexType> size_t RandomPartitionTree<IndexType>::depth() {
  return tree->depth();
}

template <typename IndexType>
IndexType RandomPartitionTree<IndexType>::find_nn(dmlc::Row<IndexType> row) { //TODO: change
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


template <typename IndexType>
void RandomPartitonTree<IndexType>::Save(dmlc::Stream *fo)
{
  tree->Save(fo);
}

template <typename IndexType>
RandomPartitionTree<IndexType>::Load(dmlc::Stream *dataFile, dmlc::Stream *treeFile)
{
  //read features
  dataFile->Read(&feature_dict);
  //read data
  data.Load(dataFile);
  //read tree
  tree = read_rptree(treeFile, &feature_dict);
}

template <typename IndexType>
RandomPartitionTree<IndexType>::Load(const char *dataFilename, const char *treeFilename)
{
  dmlc::Stream *dataFile = dmlc::Stream::Create(dataFilename, "r");
  dmlc::Stream *treeFile = dmlc::Stream::Create(treeFilename, "r");
  //read features
  dataFile->Read(&feature_dict);
  //read data
  data.Load(dataFile);
  //read tree
  tree = read_rptree(treeFile, &feature_dict);
}








};
#endif
