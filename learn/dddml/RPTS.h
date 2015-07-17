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

  real_t SparseDot(Row<IndexType> row);
};

// *** RPTSplit Definitions *** //

template <typename IndexType>
real_t RPTSplit<IndexType>::SparseDot(Row<IndexType> r1){
	size_t i,j;
	dmlc::real_t dotProduct = 0.0;
	for (i = 0, j = 0; (i < r1.length && j < feature_dict_ptr->size()); )
	{
		if (r1.index[i] == (*feature_dict_ptr)[j])
		{
			dotProduct += (r1.value[i] * d[j]);
			++i; ++j;
		}
		else if (r1.index[i] > (*feature_dict_ptr)[j])
		{
			++j;
		}
		else
		{
			++i;
		}
	}
	return dotProduct;
}

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
    //ps[i] = data[idxs[i]].SDot(d.data(), dimension);
    ps[i] = SparseDot(data[idxs[i]]);
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
  //real_t p = row.SDot(d.data(), d.size());
  real_t p = SparseDot(row);
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
  static unique_ptr<RPTNode<IndexType>> Load(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr);

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
unique_ptr<RPTNode<IndexType>> RPTNode<IndexType>::Load(dmlc::Stream *fi, std::vector<IndexType> *feature_dict_ptr)
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
    for (size_t idx : idxs) {
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

  auto find_nn(dmlc::Row<IndexType> row) -> size_t;

  auto get_rowblock() -> const dmlc::data::RowBlockContainer<IndexType> &;

  /* File IO functions */
  void Save(dmlc::Stream*);
  void Save(const char*);
  RandomPartitionTree(dmlc::Stream*, dmlc::Stream*);
  RandomPartitionTree(const char*, const char*);

private:
  std::vector<IndexType> feature_dict;
  dmlc::data::RowBlockContainer<IndexType> data;
  std::unique_ptr<rpt_impl::RPTNode<IndexType>> tree;

  inline dmlc::real_t SquareDist(const dmlc::Row<IndexType> &r1, size_t index);
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


/*
* Finds square distance between row and a datapoint index by `index'
* r1: row with global indices
* index: index to datapoint, we want to compare to
*/
template <typename IndexType>
inline dmlc::real_t RandomPartitionTree<IndexType>::SquareDist(const dmlc::Row<IndexType> &r1, size_t index)
{
	auto r2 = this->data.GetBlock()[index]; // convenient name
	size_t i,j;
	dmlc::real_t sqdist = 0.0;
	for (i = 0, j = 0; (i < r1.length && j < r2.length); )
	{
		if (r1.index[i] == feature_dict[r2.index[j]])
		{
			sqdist += (r1.value[i] - r2.value[j])*(r1.value[i] - r2.value[j]);
			++i; ++j;
		}
		else if (r1.index[i] > feature_dict[r2.index[j]])
		{
			sqdist += (r2.value[j])*(r2.value[j]);
			++j;
		}
		else
		{
			sqdist += (r1.value[i])*(r1.value[i]);
			++i;
		}
	}
	return sqdist;
}

template <typename IndexType>
size_t RandomPartitionTree<IndexType>::find_nn(dmlc::Row<IndexType> row) {
  auto leaf_idxs = tree->search(row);
  dmlc::real_t min_dist = std::numeric_limits<dmlc::real_t>::infinity();
  size_t best_idx;
  for (size_t idx : *leaf_idxs) {
  	//dmlc::real_t dist = squareDist(row, data.GetBlock()[idx]);
    dmlc::real_t dist = SquareDist(row, idx);
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
void RandomPartitionTree<IndexType>::Save(dmlc::Stream *fo)
{
  tree->Save(fo);
}

template <typename IndexType>
void RandomPartitionTree<IndexType>::Save(const char *filename)
{
  dmlc::Stream *fo = dmlc::Stream::Create(filename, "w");
  tree->Save(fo);
  delete fo;
}

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(dmlc::Stream *dataFile, dmlc::Stream *treeFile)
{
  //read features
  dataFile->Read(&feature_dict);
  //read data
  data.Load(dataFile);
  //read tree
  tree = rpt_impl::read_rptree(treeFile, &feature_dict);
}

template <typename IndexType>
RandomPartitionTree<IndexType>::RandomPartitionTree(const char *dataFilename, const char *treeFilename)
{
  dmlc::Stream *dataFile = dmlc::Stream::Create(dataFilename, "r");
  dmlc::Stream *treeFile = dmlc::Stream::Create(treeFilename, "r");
  //read features
  dataFile->Read(&feature_dict);
  //read data
  data.Load(dataFile);
  //read tree
  tree = rpt_impl::read_rptree(treeFile, &feature_dict);
  delete dataFile;
  delete treeFile;
}

};
#endif
