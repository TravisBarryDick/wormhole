#pragma once

#include "base.h"

#if DISTRIBUTED

#include "dmlc/io.h"
#include "dmlc/data.h"

#else

#include "local/io.h"
#include "local/data.h"

#endif

#include <iostream>
#include <string>
#include <cstring>
#include <random>
#include <vector>
#include <cassert>
#include <queue>
#include <memory>

namespace dddml{
using namespace dmlc;
/*
* Class for centers
*
*/

typedef std::shared_ptr<std::vector<int>> vector_int_ptr;

typedef std::shared_ptr<std::vector<std::vector<int>>> vector_vector_int_ptr;

class centers_t
{
	real_t *array;

public:
	size_t dim;
	int k;
	real_t* operator[](size_t row_id)
	{
		return &(array[row_id * dim]);
	}
	// real_t* const operator[](size_t row_id)
	// {
	// 	return &(array[row_id * dim]);
	// }
	centers_t(size_t dim, int k)
	{
		this->dim = dim;
		this-> k = k;
		array = new real_t[dim*k];
		std::memset(array, 0, sizeof(real_t)*dim*k);
	}
	~centers_t(){}
	// no automatic destructor. Call destroy() on completion
	void destroy()
	{
		if (array != NULL)
			delete[] array;
		else
			std::cout << "null pointer \n";
	}
	void reset()
	{
		std::memset(array, 0, sizeof(real_t) * dim * k);
	}
	void reset(int index, real_t val)
	{
		for (unsigned long i = index * dim; i < dim*(index + 1); ++i)
			array[i] = val;
	}

};


template<typename I>
bool update_assignments( centers_t &centers, const RowBlock<I> &data, std::vector<int> &assignments);
template<typename I>
bool update_assignments( centers_t &centers, const RowBlock<I> &data, std::vector<std::vector<int>> &assignments);
template<typename I>
void update_centers(centers_t &centers, const RowBlock<I> &data, const std::vector<int> &assignments);
template<typename I>
void update_centers(centers_t &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments);
template<typename I>
void update_one_center(centers_t &centers, const RowBlock<I> &data, const std::vector<int> &assignments, int center_id);
template<typename I>
void update_one_center(centers_t &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments, int center_id);
template<typename I>
centers_t random_init(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng);
template<typename I>
centers_t kmpp_init(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng);
template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_int_ptr assignments, centers_t &centers);
template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_vector_int_ptr assignments, centers_t &centers);
template<typename I>
std::pair<vector_int_ptr, centers_t> kmeans(const RowBlock<I> &data, int k, size_t dim, std::mt19937_64 &rng, int init = 1);
template<typename I>
std::pair<vector_vector_int_ptr, centers_t> kmeans(const RowBlock<I> &data, int k, int p, size_t dim, int center_type, std::mt19937_64 &rng, int init = 1);



template <typename I>
inline real_t squareDist(const Row<I> &r1, const Row<I> &r2)
{
	size_t i,j;
	real_t sqdist = 0.0;
	for (i = 0, j = 0; (i < r1.length && j < r2.length); )
	{
		if (r1.index[i] == r2.index[j])
		{
			sqdist += (r1.value[i] - r2.value[j])*(r1.value[i] - r2.value[j]);
			++i; ++j;
		}
		else if (r1.index[i] > r2.index[j])
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

template <typename I>
inline real_t squareDist(const Row<I> &r1, const real_t *r2, size_t dim)
{
	#if 0
	real_t sqdist = 0;
	for (int i = 0; i < this->dim; ++i) sqdist += r2[i] * (r2[i]);
	for (int i = 0; i < r1.length; ++i) sqdist += (r1.value[i])*(r1.value[i]);
	sqdist -= 2 * r1.SDot(r2, this->dim);
	return sqdist;
	#endif
	//#if 0
	CHECK(r2 != NULL);
	size_t i,j;
	real_t sqdist = 0.;
	for (i = 0, j = 0; (i < r1.length && j < dim); )
	{
		if (r1.index[i] == j)
		{
			sqdist += (r1.value[i] - r2[j])*(r1.value[i] - r2[j]);
			++i; ++j;
		}
		else if (r1.index[i] > j)
		{
			sqdist += (r2[j] * r2[j]);
			++j;
		}
		else
		{
			CHECK(-2 == 0);
			//sqdist += (r1.value[i] * r1.value[i]);
			//++i;
		}
	}
	return sqdist;
	//#endif
}

inline real_t squareDist(const real_t *array1, const real_t *array2, size_t dim)
{
	CHECK(array1 != NULL);
	CHECK(array2 != NULL);
	real_t sqdist  = 0;
	for (int i = 0; i < dim; ++i)
	{
		sqdist += (array1[i] - array2[i]) * (array1[i] - array2[i]);
	}
	return sqdist;
}

inline void reset(real_t *array, size_t dim)
{
	CHECK(array != NULL);
	for (int i = 0; i < dim; ++i) array[i] = 0.0;
}

template <typename I>
inline void add_into(real_t *arr, size_t dim, const Row<I> &r1)
{
	CHECK(arr != NULL);
	CHECK(r1.index[r1.length - 1]	 < dim);
	for(size_t i = 0; i < r1.length; ++i)
	{
		arr[r1.index[i]] += r1.weight * r1.value[i];

	}
}

/* divide *this by a scalar div */
inline void divide_by(real_t *arr, size_t dim, const real_t div)
{
	CHECK(arr != NULL);
	for(size_t i = 0; i < dim; ++i)
	{
		arr[i] = (div == 0) ? 0: arr[i] / div;
	}
}

inline real_t *initialize_center(size_t dim)
{
	real_t *array = new real_t[dim];
	for (int i = 0; i < dim; ++i) array[i] = 0.0;
	return array;
}

inline real_t **initialize_centers(size_t dim, int k)
{
	real_t **centers = new real_t*[k];
	for (int i = 0; i < k; ++i)
	{
		centers[i] = new real_t[dim];
		std::memset(centers[i], 0, sizeof(real_t)*dim);
	}
	return centers;
}

// inline void destroy_center(real_t *arr)
// {
// 	CHECK(arr != NULL);
// 	delete[] arr;
// }

// inline void destroy_centers(real_t **centers, int k)
// {
// 	CHECK(centers != NULL);
// 	for (int i = 0; i < k; ++i)
// 	{
// 		CHECK(centers[i] != NULL);
// 		delete[] centers[i];
// 	}
// 	delete[] centers;
// }


/*
* functions to find nearest point from a given collection
*/

template<typename I>
int find_closest(const Row<I> &row, centers_t &centers)
{
	size_t dim = centers.dim; int k = centers.k;
	if (k == 0) return -1;
	else if (k == 1) return 0;
	int min_index = 0;
	//std::cout << squareDist(row, centers[0], dim) << std::endl;
	real_t min_dist = squareDist(row, centers[0], dim);
	real_t cur_dist;
	for (size_t i = 1; i < k; ++i)
	{
		cur_dist = squareDist(row, centers[i], dim);
		if (cur_dist < min_dist)
		{
			min_dist = cur_dist;
			min_index = i;
		}
	}
	return min_index;
}


/*
* Find closest point from vector ignoring index i
*/
int find_closest(centers_t &all_centers, int center)
{
	size_t dim = all_centers.dim; int k = all_centers.k;
	int min_index = (center != 0) ? 0 : 1;
	real_t min_dist = squareDist(all_centers[min_index], all_centers[center], dim);
	real_t cur_dist = 0.;
	for (size_t i = 1 + min_index; i < k; ++i)
	{
		if (i == center)
		{
			continue;
		}
		else
		{
			cur_dist = squareDist(all_centers[i], all_centers[center], dim);
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist;
				min_index = i;
			}
		}
	}
	return min_index;
}

// returns p nearest centers from nearest to farthest
template<typename I>
int *find_p_closest(int p, const Row<I> &row,  centers_t &centers)
{
	size_t dim = centers.dim; int k = centers.k;
	if (k < p)
	{
	//should not occur:
		exit(0);
	}
	else
	{
		real_t sqdist;
		//heap select to efficiently find p things
		std::priority_queue<std::pair<real_t, int>> pq;
		for (int i = 0; i < p; ++i) //insert for p elements into the heap
		{
			pq.push(std::pair<real_t, int>(squareDist(row, centers[i], dim), i));
		}
		//insert rest of the distances while maintaing p smallest in the heap
		for (int i = p; i < k; ++i)
		{
			sqdist = squareDist(row, centers[i], dim);
			if (pq.top().first > sqdist)
			{
				pq.pop();
				pq.push(std::pair<real_t, int>(sqdist, i));
			}
		}
		//read off assignments in reverse order
		int *assignments = new int[p];
		std::memset(assignments, 0, sizeof(int) * p);
		for (int i = p-1; i >= 0; --i)
		{
			assignments[i] = pq.top().second;
			pq.pop();
		}
		return assignments;
	}
}

///////////////////////////////////////////////
// Saving and Reading Clustering Assignments //
///////////////////////////////////////////////

#if DISTRIBUTED

void save_assignments(dmlc::Stream *fo,
                      const std::vector<std::vector<int>> *assignments) {
  size_t n = assignments->size();
  fo->Write(&n, sizeof(size_t));
  for (size_t i = 0; i < n; ++i) {
    fo->Write((*assignments)[i]);
  }
}

void save_assignments(const char *outpath,
                      const std::vector<std::vector<int>> *assignments) {
  auto outstream = Stream::Create(outpath, "w");
  save_assignments(outstream, assignments);
  delete outstream;
}

void read_assignments(dmlc::Stream *fi,
                      std::vector<std::vector<int>> *assignments) {
  size_t n;
  fi->Read(&n, sizeof(size_t));
  assignments->resize(n);
  for (size_t i = 0; i < n; ++i) {
    fi->Read(&(*assignments)[i]);
  }
}

void read_assignments(const char *inpath,
                      std::vector<std::vector<int>> *assignments) {
  auto instream = Stream::Create(inpath, "r");
  read_assignments(instream, assignments);
  delete instream;
}
#endif

}//namespace dddml
