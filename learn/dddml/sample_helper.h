#pragma once
#include <string>
#include <cstring>
#include <random>
namespace dddml{

// HELPER CLASSES/FUNCTIONS
inline void swap__(int *a, int pos1, int pos2)
{
	int temp = a[pos1];
	a[pos1] = a[pos2];
	a[pos2] = temp;
}
/*
* Sample "size" number of elements from 1:range
*/
int *SampleWithoutReplacement(unsigned size, unsigned range, std::mt19937_64 &rng)
{
	
	int arr [range];
	int *ret = new int[size];
	for (unsigned i = 1; i <= range; ++i) arr[i] = i;
	
	int last = range; 
	int rnd;
	std::uniform_int_distribution<int> dunif(0, last-1);
	for (unsigned i = 0; i < size; ++i)
	{
		rnd = dunif(rng);
		ret[i] = arr[rnd];
		swap__(arr, rnd, last - 1);
		--last;
	}
	return ret;
}

/*
* Sample i from 0 to (l-1) w.p. proportional to weights[i]
*/
inline int weightedSample(float *weights, size_t length, std::mt19937_64 &rng)
{
	std::discrete_distribution<> gen(weights, weights+length);
	return gen(rng);
}

}//namespace dddml
