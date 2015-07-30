#include <string>
#include <random>
#include <vector>
#include "kmeans_helper.h"
#include "sample_helper.h"
#include <memory>
#include <queue>
#include <array>
#include <algorithm>
#include <sstream>
#include "RPTS.h"

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
*	IMPLEMENTATION of k-MEANS with k-means++ INITIALIZATION
* -----------------------------------------------------------
*/

namespace dddml{
using namespace dmlc;
using namespace dmlc::data;




/*
* Update assignment step of k-means based on Voronoi partitions
* /param centers: vector of centers
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
* /return true if any assignments have changed. False otherwise
*/
template<typename I>
bool update_assignments(centers_t &centers, const RowBlock<I> &data, std::vector<int> &assignments)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t n = data.size;
	CHECK(data.size == assignments.size());
	bool changed = false;
	int assignment;
	for (size_t i = 0; i < n; ++i)
	{

		assignment = find_closest(data[i], centers);//, dim, k);
		if (assignment != assignments[i])
		{
			assignments[i] = assignment;
			changed = true;
		}
	}
	return changed;
}

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
* /param centers: vector of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
*/

template<typename I>
void update_centers(centers_t &centers, const RowBlock<I> &data, const std::vector<int> &assignments)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t n = data.size;
	CHECK(data.size == assignments.size());
	centers.reset();
	int counts[k];
	for (int i = 0; i < k; ++i) counts[i] = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		if (assignments[i] >= 0)
		{ //valid assignment
			add_into(centers[assignments[i]], dim, data[i]);
			counts[assignments[i]] += (data.weight == NULL) ? 1 : data.weight[i];
		}
	}
	//average
	for (int i = 0; i < k; ++i)
	{
		divide_by(centers[i], dim, counts[i]);
	}
}


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
* /param centers: vector of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
* /param center_id: ID of center to recompute
*/

template<typename I>
void update_one_center(centers_t &centers, const RowBlock<I> &data, const std::vector<int> &assignments, int center_id)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t n = data.size;
	CHECK(data.size == assignments.size());

	reset(centers[center_id], dim);

	int count = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		if (assignments[i] == center_id)
		{
			add_into(centers[center_id], dim, data[i]);
			count += (data.weight == NULL) ? 1 : data.weight[i];
		}
	}
	//average
	divide_by(centers[center_id], dim, count);
}


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


/*
* Pick up random centers to initialize k-means (code: 0)
*/
template<typename I>
centers_t random_init(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng)
{
	size_t numData = data.size;
	int *sample = SampleWithoutReplacement(k, numData, rng);
	//alternative 1:
	centers_t centers(dim, k);
	for (int i = 0; i < k; ++i)
	{
		add_into(centers[i], dim, data[sample[i]]);
	}

	return centers;
}

/*
*	k-means++ initialization (code: 1)
*/
template<typename I>
centers_t kmpp_init(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng)
{
	real_t weight;
	centers_t centers(dim, k);

	size_t numData = data.size;
	real_t sqdists[numData]; //distances
	for (int i = 0; i < numData; ++i) sqdists[i] = 0;
	// initialize first center
	std::uniform_int_distribution<> dis(0, numData-1);

	add_into(centers[0], dim, data[dis(rng)]);

	//initialize distances
	for (int j = 0; j < numData; ++j)
	{
		weight = (data.weight == NULL) ? 1 : data.weight[j];
		sqdists[j] = (weight) * squareDist(data[j], centers[0], dim);
	}
	//loop for the next (k-1) centers
	for (int i = 1; i < k; ++i)
	{
		//sample next center
		int next_center = weightedSample(sqdists, numData, rng);
		add_into(centers[i], dim, data[next_center]);

		//update distances
		for (int j = 0; j < numData; ++j)
		{
			weight = (data.weight == NULL) ? 1 : data.weight[j];
			sqdists[j] = std::min(sqdists[j], squareDist(data[j], centers[i], dim));
		}
	}
	return centers;
}

template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_int_ptr assignments, centers_t &centers)
{
	size_t dim = centers.dim; int k = centers.k;
	size_t n = data.size;
	real_t obj = 0;
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		if ((*assignments)[i] != -1)
		{
			++ count;
			real_t weight = (data.weight == NULL) ? 1 : data.weight[i];
			obj += weight * squareDist(data[i], centers[(*assignments)[i]], dim);
		}
	}
	if (count != n) std::cerr << "ERROR@!\n";
	return obj/count;
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
std::pair<vector_int_ptr, centers_t> kmeans(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng, int init )
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
	vector_int_ptr assignments = std::make_shared<std::vector<int>>(data.size, 0);
	start = GetTime();
	int iter_count = 0;
	while(changed)
	{
		double start1 = GetTime();
		changed = update_assignments(centers, data, *assignments);//, dim, k);
		//std::cout << "\tupdate Assignment: " << GetTime() - start1 << std::endl;

		int counts[k] ;
		for (int i = 0; i < k; ++i) counts[i] = 0;
		for (int i = 0; i < assignments->size(); ++i)
		{
			++counts[(*assignments)[i]];
		}
		std::cout << '\t';
		for (int i = 0; i < k; ++i) std::cout << counts[i] << ' ';
		std::cout << std::endl;

		start1 = GetTime();
		update_centers(centers, data, *assignments);
		//std::cout << "\tupdate centers: " << GetTime() - start1 << std::endl;
		++iter_count;
		std::cout << iter_count << ": objective: " << kmeans_objective(data, assignments, centers) << std::endl;
	}
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Finished k-means: " << iter_count << " iterations in " << time << " sec";
	#else
	std::cout << "Finished k-means: " << iter_count << " iterations in " << time << " sec" << std::endl;
	#endif

	return std::pair<vector_int_ptr, centers_t>(assignments, centers);
}

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


#if !DISTRIBUTED
template<typename I>
void save_data_to_file(const char * filename, const RowBlock<I> &data, vector_int_ptr assignments, Stream *fo = NULL)
{
	libsvmwrite<I>(filename, data, *assignments);
}

template<typename I>
void save_data_to_file(const char *filename, const RowBlock<I> &data, vector_vector_int_ptr assignments, Stream *fo = NULL)
{
	libsvmwrite(filename, data, *assignments);
}


template<typename I>
void save_centers_to_file(Stream *fo, centers_t &centers)
{
	//TODO
}
#endif

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

template <typename I>
int merge_and_split(const RowBlock<I> &data, vector_int_ptr assignments, centers_t &centers, real_t lower_bound, real_t upper_bound, int k, size_t dim, std::mt19937_64 &rng)
{
	std::stringstream outs;
	//preprocess
	int current_k = k;
	real_t current_lower = lower_bound,
			current_upper = upper_bound;
	int counts[k] ;
	for (int i = 0; i < k; ++i) counts[i] = 0;
	int nmerge = 0,  nsplit = 0;
	for (int i = 0; i < assignments->size(); ++i)
	{
		++counts[(*assignments)[i]];
	}
	std::vector<int> permutation(k);
	std::vector<int> deleted_indices;
	outs << current_k << std::endl;
	#if DISTRIBUTED
	LOG(INFO) << outs;
	#else
	std::cout << outs << std::endl;
	#endif
	outs.str(std::string());
	//merge
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] < current_lower)
		{
			//cluster too small. Need to merge it with nearest cluster.
			deleted_indices.push_back(i);
			int new_index = find_closest(centers, i);
			++nmerge;
			counts[new_index] += counts[i];
			counts[i] = -1; //destroyed cluster
			--current_k;
			outs << "Merged cluster " << i << " with cluster " << new_index << std::endl;
			#if DISTRIBUTED
			LOG(INFO) << outs;
			#else
			std::cout << outs << std::endl;
			#endif
			outs.str(std::string());
			for (int j = 0; j < assignments->size(); ++j)
			{
				if ((*assignments)[j] == i) (*assignments)[j] = new_index;
			}
			//update center:
			update_one_center(centers, data, (*assignments), new_index);
		}
	}
	//modify upper bound (it gets looser) and leave lower bound intact
	current_upper *= k / current_k;



	//split
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	int current_index = k; //starting index for new clusters formed
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] > current_upper)
		{
			//cluster too big. Need to split.
			nsplit++;
			int num_new_clusters_for_this_cluster = counts[i] / current_upper + (counts[i] != current_upper); //ceil
			outs << "split cluster " << i << " to clusters " << current_index  << " ... " << (current_index + num_new_clusters_for_this_cluster - 1) << std::endl;
			#if DISTRIBUTED
			LOG(INFO) << outs;
			#else
			std::cout << outs << std::endl;
			#endif
			outs.str(std::string());
			deleted_indices.push_back(i); //delete this index
			std::vector<int> ones = _ones(num_new_clusters_for_this_cluster);
			std::discrete_distribution<int> dis (ones.begin(), ones.end());
			for (int j = 0; j < assignments->size(); ++j)
			{
				if ((*assignments)[j] == i)  //ramdomly assign to one of the split clusters
				{
					int temp = dis(rng);
					(*assignments)[j] = temp + current_index;
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
		int cl = (*assignments)[i];
		++counts1[cl];
	}
	for (int i = 0; i < current_index; ++i) {std::cout << i << ": " << counts1[i] << std::endl; }


	//clean-up
	//re-map indices >=current_k to deleted indices
	for (int j = 0; j < assignments->size(); ++j)
	{
		int temp = (*assignments)[j] - current_k;
		if (temp >= 0) (*assignments)[j] = deleted_indices[temp];
	}
	#if DISTRIBUTED
	LOG(INFO) << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k;
	#else
	std::cout << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k << std::endl;
	#endif
	return current_k;
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
	std::cout << current_k << std::endl;
	//merge
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] < current_lower)
		{
			//cluster too small. Need to merge it with nearest cluster.
			deleted_indices.push_back(i);
			int new_index = find_closest(centers, i);
			++nmerge;
			counts[new_index] += counts[i];
			counts[i] = -1; //destroyed cluster
			--current_k;
			std::cout << "Merged cluster " << i << " with cluster " << new_index << std::endl;
			for (int j = 0; j < assignments->size(); ++j)
			{
				for (int h = 0; h < p; ++h)
					if ((*assignments)[j][h] == i) (*assignments)[j][h] = new_index;
			}
			//update center:
			update_one_center(centers, data, (*assignments), new_index);
		}
	}
	//modify upper bound (it gets looser) and leave lower bound intact
	current_upper *= k / current_k;



	//split
	for (int i = 0; i < k; ++i) permutation[i] = i;
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
						//std::cout << '*' << temp << std::endl;
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
	for (int i = 0; i < current_index; ++i) {std::cout << i << ": " << counts1[i] << std::endl; }


	//clean-up
	//re-map indices >=current_k to deleted indices
	for (int j = 0; j < assignments->size(); ++j)
	{
		for (int h = 0; h < p; ++h)
		{
			int temp = (*assignments)[j][h] - current_k;
			if (temp >= 0) (*assignments)[j][h] = deleted_indices[temp];
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


int main(int argc, char *argv[])
{
	using namespace dddml;
	using namespace std;
	using namespace dmlc;
	ArgParser parser;
	if (argc > 1 && strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
	parser.ReadArgs(argc - 2, argv + 2);
	dddml::dddmlConfig conf;
	parser.ParseToProto(&conf);


	std::random_device rd;
	std::mt19937_64 rng(rd());

	/*
		Parameters:
			- n_clusters: K: number of clusters
			- replication_factor: p: replication factor
			- cluster_lower_bound: l: lower bound (fraction)
			- cluster_upper_bound: L: upper bound (fraction)
			- data_file: file with idx_map and localized row block
			- assignments_file: file to write assignments to
	*/

	int k = conf.n_clusters(),
		p = conf.replication_factor();
	real_t lfrac = conf.cluster_lower_bound(),
		Lfrac = conf.cluster_upper_bound();
	const char *data_file = conf.sample_filename().c_str(), *out_file = conf.assignments_filename().c_str();

	auto readpair = readSamplingOutput(data_file);
	auto idx_dict = readpair.first;
	int dim = idx_dict->size(); // dim = size of idx_dict
	dmlc::data::RowBlockContainer<FeaID>* data_rbc = readpair.second;  // pointer to row block container with data
	auto data = data_rbc->GetBlock();
	int n = data.size;

	real_t lb = lfrac * n * p / k,
			ub = Lfrac * n * p / k;

	std::cout << "DEBUG: data_rbc->Size() = " << data_rbc->Size() << std::endl;
	for (size_t i = 0; i < data_rbc->Size(); ++i) {
		std::cout << data[i].length << " ";
	}
	std::cout << std::endl;

	auto output = kmeans(data,  k, p , dim, rng, 0);
	auto assignments = output.first;
	auto centers = output.second;

	//printing
	int counts[k];
	for (int i = 0; i < k; ++i) {counts[i] = 0;}

	for (int i = 0; i < assignments->size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			int cl = (*assignments)[i][j];
			++counts[cl];
		}
	}
	for (int i = 0; i < k; ++i)
	{
		std::cout << i << ": " << counts[i] << "; " << std::endl;
	}
	std::cout << "-----------------------\n";


	//heuristics
	int new_k = merge_and_split(data, assignments, centers, lb, ub, k, dim, rng);

	//printing
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
	centers.destroy();

	// Save the number of clusters and the assignments file
	std::ofstream num_clusters(conf.num_clusters_filename().c_str(), std::ios::out);
	num_clusters << new_k << std::endl;
	num_clusters.close();
	save_assignments(out_file, &(*assignments));
	
	LOG(INFO) << "FINISHED CLUSTERING";
	
	// Build a random partition tree on the sample and save it to file
	RandomPartitionTree<FeaID> rpt(rng, static_cast<int>(conf.dimension()),
																 static_cast<int>(conf.n_0()), *data_rbc,
																 *idx_dict);
	rpt.Save(conf.rpt_filename().c_str());
}

#else

int main1()
{
	using namespace dddml;
	using namespace std;
	std::random_device rd;
	std::mt19937_64 rng(rd());
	dmlc::data::RowBlockContainer<int> rbc = dmlc::data::libsvmread("./mnist.txt");
	RowBlock<int> block = rbc.GetBlock();
	int n = block.size;
	auto rb = block;
	int k = 10;
	size_t dim = 784;


	//#if 0
	auto output = kmeans(block,  k, /*p */ /*dim */ dim, rng, 1);
	std::cout << "-----------------------\n";
	vector_int_ptr assignments = output.first;
	auto centers = output.second;
	//for (auto i : *assignments)
	//	std::cout << i << std::endl;
	int counts[k];
	real_t distances[k];
	for (int i = 0; i < k; ++i) {counts[i] = 0; distances[i] = 0;}

	for (int i = 0; i < assignments->size(); ++i)
	{
		int cl = (*assignments)[i];
		++counts[cl];
		distances[cl] += squareDist(block[i], centers[cl], dim);
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
		int cl = (*assignments)[i];
		++counts1[cl];
	}
	for (int i = 0; i < new_k; ++i)
	{
		std::cout << i << ": " << counts1[i] << "; " << std::endl;
	}
	//#endif
	centers.destroy();
	save_data_to_file("./sample_out", block, assignments);

}

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
