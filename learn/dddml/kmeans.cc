#include <string>
#include <random>
#include <vector>
#include "kmeans.h"
#include <memory>
#include <queue>
#include <array>
#include <algorithm>
#include <limits>
#include <sstream>
#include "RPTS.h"
#include "config_tools.h"

#if DISTRIBUTED

#include <dmlc/data.h>
#include <dmlc/logging.h>
#include <dmlc/io.h>
#include "data/row_block.h"
#include <dmlc/timer.h>
#include "dddml_config.pb.h"
#include "base/arg_parser.h"
#include "ps.h"

#else

#include <iostream>
#include "local/data.h"
#include "local/io.h"
#include "local/row_block.h"
#include "local/timer.h"
#include "local/libsvm_parser.h"

#endif



/*
* -----------------------------------------------------------
*	IMPLEMENTATION of fault tolerant k-MEANS with k-means++ INITIALIZATION
* -----------------------------------------------------------
*/

namespace dddml{
using namespace dmlc;
using namespace dmlc::data;




/*
* Update assignment step of k-means based on Voronoi partitions
* /param centers: data structure of centers
* /param data: datapoints that are to be reassigned
* /param assignments: vector of p-vector of assignments from 1 to k. It will be modified
* /return true if any assignments have changed. False otherwise
*/


template<typename I>
bool update_assignments(centers_t &centers, const RowBlock<I> &data, std::vector<std::vector<int>> &assignments)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t  p = assignments[0].size(), n = data.size;
	CHECK(data.size == assignments.size());
	bool changed = false;
	std::vector<int> assignment;
	for (size_t i = 0; i < n; ++i)
	{
		auto p_closest = find_p_closest(p, data[i], centers);
		for (size_t j = 0; j < p; ++j)
		{
			if (assignments[i][j] != p_closest[j]){
				changed = true;
				assignments[i][j] = p_closest[j];
			}
		}
		delete[] p_closest;
	}
	return changed;
}

/*
* Update centers step of k-means
* /param centers: data structure of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of p-vector of assignments from 1 to k. It will be modified
*/
template<typename I>
void update_centers(centers_t &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t p = assignments[0].size(), n = data.size;
	CHECK(data.size == assignments.size());

	centers.reset();

	int counts[k];
	for (int i = 0; i < k; ++i) counts[i] = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if (assignments[i][j] >= 0)
			{ //valid assignment
				add_into(centers[assignments[i][j]], dim, data[i]);
				counts[assignments[i][j]] += (data.weight == NULL) ? 1 : data.weight[i];
			}
		}
	}
	//average
	for (int i = 0; i < k; ++i)
	{
		divide_by(centers[i], dim, counts[i]);
		//centers[i]->updateNorm();
	}
}

/*
* Update a single center: for balancing heuristics
* /param centers: data structure of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of p-vector assignments from 1 to k. It will be modified
* /param center_id: ID of center to recompute
*/



template<typename I>
void update_one_center(centers_t &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments, int center_id)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t p = assignments[0].size(), n = data.size;
	CHECK(data.size == assignments.size());
	reset(centers[center_id], dim);
	int count = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if (assignments[i][j] == center_id)
			{ //valid assignment
				add_into(centers[center_id], dim, data[i]);
				count += (data.weight == NULL) ? 1 : data.weight[i];
			}
		}
	}
	//average
	divide_by(centers[center_id], dim, count);
}


template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_vector_int_ptr assignments, centers_t &centers)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t n = data.size;
	int p = (*assignments)[0].size();
	int count = 0;
	real_t obj = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if ((*assignments)[i][j] != -1)
			{
				++count;
				real_t weight = (data.weight == NULL) ? 1 : data.weight[i];
				obj += weight * squareDist(data[i], centers[(*assignments)[i][j]], dim);
			}
		}
	}
	return obj / count;
}


/*
* -----------------------------------------------------------------------
*	k-means
* ----------------------------------------------------------------------
*/

template<typename I>
std::pair<vector_vector_int_ptr, centers_t> kmeans(const RowBlock<I> &data, int k, int p, size_t dim, std::mt19937_64 &rng, int init)
{
	double start, time;

	#if DISTRIBUTED
	LOG(INFO) << "Starting k-means";
	#else
	std::cout << "Starting k-means" << std::endl;
	#endif
	#if DISTRIBUTED
	LOG(INFO) << "Initializing centers..";
	#else
	std::cout << "Initializing centers.." << std::endl;
	#endif

	centers_t centers(dim, k);
	start = GetTime();
	if (init == 0)
		centers = random_init(data, k, dim, rng);
	else //kmpp
		centers = kmpp_init(data, k, dim, rng);
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Initialization done in " << time << " sec";
	#else
	std::cout << "Initialization done in " << time << " sec" << std::endl;
	#endif


	#if DISTRIBUTED
	LOG(INFO) << "Starting Lloyd's iterations";
	#else
	std::cout << "Starting Lloyd's iterations" << std::endl;
	#endif

	bool changed = true;
	vector_vector_int_ptr assignments = std::make_shared<std::vector<std::vector<int>>>(data.size, std::vector<int>(p));
	start = GetTime();
	int iter_count = 0;
	while(changed)
	{
		changed = update_assignments(centers, data, *assignments);

		int counts[k] ;
		for (int i = 0; i < k; ++i) counts[i] = 0;
		for (int i = 0; i < assignments->size(); ++i)
		{
			for (int j = 0; j < p; ++j)
				++counts[(*assignments)[i][j]];
		}
		std::cout << '\t';
		for (int i = 0; i < k; ++i) std::cout << counts[i] << ' ';
		std::cout << std::endl;

		update_centers(centers, data, *assignments);
		++iter_count;
		std::cout << iter_count << ": objective: " << kmeans_objective(data, assignments, centers) << std::endl;
	}
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Finished k-means: " << iter_count << " iterations in " << time << " sec";
	#else
	std::cout << "Finished k-means: " << iter_count << " iterations in " << time << " sec" << std::endl;
	#endif

	return std::pair<vector_vector_int_ptr, centers_t>(assignments, centers);
}

/*
*	MERGE AND SPLIT HEURISTICS:
*	 Heuristics to merge clusters that are too small and split clusters that are too big
*
*/

std::vector<int> _ones(int len)
{
	std::vector<int> ones(len);
	for (int i = 0; i < len; ++i) ones[i] = 1;
	return ones;
}


/* p > 1 */
template <typename I>
int merge_and_split(const RowBlock<I> &data, vector_vector_int_ptr assignments, centers_t &centers, real_t lower_bound, real_t upper_bound, int k, size_t dim, std::mt19937_64 &rng)
{
	//preprocess
	int current_k = k, p = (*assignments)[0].size();
	real_t current_lower = lower_bound,
			current_upper = upper_bound;
	int counts[k] ;
	for (int i = 0; i < k; ++i) counts[i] = 0;
	int nmerge = 0,  nsplit = 0;
	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		++counts[(*assignments)[i][j]];
	}
	std::vector<int> permutation(k);
	std::vector<int> deleted_indices;
	//merge
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] < current_lower)
		{
			//cluster i is too small. Need to merge it with nearest cluster.
			deleted_indices.push_back(i);
			int new_index = find_closest(centers, i);
			++nmerge;
			counts[new_index] += counts[i];
			counts[i] = -1; //destroyed cluster
			--current_k;
			std::cout << "Merged cluster " << i << " with cluster " << new_index << std::endl;
			int counter = 0;
			for (int j = 0; j < assignments->size(); ++j)
			{
				for (int h = 0; h < p; ++h){
					if ((*assignments)[j][h] == i) {(*assignments)[j][h] = new_index;}
				}
			}
			//update centers:
			centers.reset(i, std::numeric_limits<real_t>::max());
			update_one_center(centers, data, (*assignments), new_index);
		}
	}
	//modify upper bound (it gets looser) and leave lower bound intact
	current_upper *= k / current_k;



	//split
	for (int i = 0; i < k; ++i) permutation[i] = i; //permute
	std::shuffle(permutation.begin(), permutation.end(), rng);
	int current_index = k; //starting index for new clusters formed
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] > current_upper)
		{
			//cluster too big. Need to split.
			nsplit++;
			int num_new_clusters_for_this_cluster = counts[i] / current_upper + (counts[i] != current_upper); //ceil
			std::cout << "split cluster " << i << " to clusters " << current_index  << " ... " << (current_index + num_new_clusters_for_this_cluster - 1) << std::endl;
			deleted_indices.push_back(i); //delete this index
			std::vector<int> ones = _ones(num_new_clusters_for_this_cluster);
			std::discrete_distribution<int> dis (ones.begin(), ones.end());
			for (int j = 0; j < assignments->size(); ++j)
			{
				for (int h = 0; h < p; ++h)
				{
					if ((*assignments)[j][h] == i)  //ramdomly assign to one of the split clusters
					{
						int temp = dis(rng);
						(*assignments)[j][h] = temp + current_index;
					}
				}
			}
			current_index += num_new_clusters_for_this_cluster; //next cluster starts from this index
			/*
			if we split into t clusters, t-1 new clusters are formed
			*/
			current_k += (num_new_clusters_for_this_cluster - 1);
		}
	}

	//pre-clean-up
	int counts1[current_index];
	for (int i = 0; i < current_index; ++i) {counts1[i] = 0; }

	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			int cl = (*assignments)[i][j];
			++counts1[cl];
		}
	}
	//for (int i = 0; i < current_index; ++i) {std::cout << i << ": " << counts1[i] << std::endl; }


	//clean-up
	int swap_map[current_index];
	bool to_be_swapped[current_index];
	for (int i = 0; i < current_index; ++i) {to_be_swapped[i] = 0; swap_map[i] = -1;}
	for (int ii = 0, jj = current_index - 1; ii < jj ; )
	{
		while(counts1[ii] > 0) ++ii;
		while(counts1[jj] <= 0) --jj;
		if (ii >= jj ) break;
		swap_map[jj] = ii;
		to_be_swapped[jj] = true;
		++ii; --jj;
	}
	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			int temp = (*assignments)[i][j];
			if (to_be_swapped[temp]) (*assignments)[i][j] = swap_map[temp];
		}
	}

	#if DISTRIBUTED
	LOG(INFO) << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k;
	#else
	std::cout << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k << std::endl;
	#endif
	return current_k;
}



} //namespace dddml

#if DISTRIBUTED
namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps


using FeaID = unsigned long long;
using myPair = std::pair<std::vector<FeaID>*, dmlc::data::RowBlockContainer<FeaID>*>;

myPair readSamplingOutput(const char *filename)
{
	dmlc::Stream *output = dmlc::Stream::Create(filename, "r");
	std::vector<FeaID> *idx = new std::vector<FeaID>();
	dmlc::data::RowBlockContainer<FeaID> *data = new dmlc::data::RowBlockContainer<FeaID>();
	// read indices
	bool read = output->Read(idx);
	// read rowblock
	data->Load(output);
	return myPair(idx, data);
}

int main(int argc, char *argv[]) {
  using namespace dddml;
  using namespace std;
  using namespace dmlc;

  ConfigWrapper cfg(argv[1]);

  std::mt19937_64 rng(cfg.clustering_seed);
  int k = cfg.clustering_num_clusters;
  int p = cfg.clustering_replication;
  real_t lfrac = cfg.clustering_lower_capacity;
  real_t Lfrac = cfg.clustering_upper_capacity;
  const char *data_file = cfg.dispatch_sample_file.c_str();
  const char *out_file = cfg.dispatch_assignments_file.c_str();

  auto readpair = readSamplingOutput(cfg.dispatch_sample_file.c_str());
  auto idx_dict = readpair.first;
  int dim = idx_dict->size();  // dim = size of idx_dict
  dmlc::data::RowBlockContainer<FeaID> *data_rbc =
      readpair.second;  // pointer to row block container with data
  auto data = data_rbc->GetBlock();
  int n = data.size;

  real_t lb = lfrac * n * p / k;
  real_t ub = Lfrac * n * p / k;

  auto output = kmeans(data, k, p, dim, rng, 0);
  auto assignments = output.first;
  auto centers = output.second;

  // printing
  int counts[k];
  for (int i = 0; i < k; ++i) {
    counts[i] = 0;
  }

  for (int i = 0; i < assignments->size(); ++i) {
    for (int j = 0; j < p; ++j) {
      int cl = (*assignments)[i][j];
      ++counts[cl];
    }
  }
  for (int i = 0; i < k; ++i) {
    std::cout << i << ": " << counts[i] << "; " << std::endl;
  }
  std::cout << "-----------------------\n";

  // heuristics
  int new_k = merge_and_split(data, assignments, centers, lb, ub, k, dim, rng);

  // printing
  int counts1[new_k];
  for (int i = 0; i < new_k; ++i) {
    counts1[i] = 0;
  }

  for (int i = 0; i < assignments->size(); ++i) {
    for (int j = 0; j < p; ++j) {
      int cl = (*assignments)[i][j];
      ++counts1[cl];
    }
  }
  for (int i = 0; i < new_k; ++i) {
    std::cout << i << ": " << counts1[i] << std::endl;
  }
  centers.destroy();

  // Save the number of clusters and the assignments file
  std::ofstream num_clusters(cfg.dispatched_num_clusters_file.c_str(),
                             std::ios::out);
  num_clusters << new_k << std::endl;
  num_clusters.close();
  save_assignments(out_file, &(*assignments));

  LOG(INFO) << "FINISHED CLUSTERING";

  // Build a random partition tree on the sample and save it to file
  RandomPartitionTree<FeaID> rpt(rng, static_cast<int>(dim),
                                 cfg.dispatch_rpt_n0, *data_rbc, *idx_dict);
  rpt.Save(cfg.dispatch_rpt_file.c_str());
}

#else

int main()
{
	using namespace dddml;
	using namespace std;
	std::random_device rd;
	std::mt19937_64 rng(rd());
	dmlc::data::RowBlockContainer<int> rbc = dmlc::data::libsvmread("./mnist.txt");
	RowBlock<int> block = rbc.GetBlock();
	int n = block.size;
	auto rb = block;
	int k = 10, p = 1;
	size_t dim = 784;


	//#if 0
	auto output = kmeans(block,  k, /*p */ p , /*dim */ dim, rng, 1);
	std::cout << "-----------------------\n";
	auto assignments = output.first;
	auto centers = output.second;
	//for (auto i : *assignments)
	//	std::cout << i << std::endl;
	int counts[k];
	real_t distances[k];
	for (int i = 0; i < k; ++i) {counts[i] = 0; distances[i] = 0;}

	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			int cl = (*assignments)[i][j];
			++counts[cl];
		}
		//distances[cl] += squareDist(block[i], centers[cl], dim);
	}
	for (int i = 0; i < k; ++i)
	{
		std::cout << i << ": " << counts[i] << "; " << std::endl;
	}

	std::cout << "-----------------------\n";
	real_t lb = 0.5 * n / k , ub = 2 * n / k;
	std::cout << "Bounds: " << lb << ' ' << ub << std::endl;

	int new_k = merge_and_split(block, assignments, centers, lb, ub, k, dim, rng);

	int counts1[new_k];
	for (int i = 0; i < new_k; ++i) {counts1[i] = 0; }

	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			int cl = (*assignments)[i][j];
			++counts1[cl];
		}
	}
	for (int i = 0; i < new_k; ++i)
	{
		std::cout << i << ": " << counts1[i] << "; " << std::endl;
	}
	//#endif
	centers.destroy();
	save_data_to_file("./sample_out", block, assignments);

}


#endif
