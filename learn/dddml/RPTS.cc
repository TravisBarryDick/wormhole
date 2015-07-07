#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <sstream>

#include "RPTS.h"
#include "local/io.h"
#include "local/libsvm_parser.h"
#include "kmeans_helper.h"

using namespace std;
using namespace std::placeholders;
using namespace dddml;
using namespace dddml::rpt_impl;
using namespace dmlc;
using namespace dmlc::data;

auto query_dist(function<int(Row<int>)> f, const RowBlockContainer<int> &data,
                Row<int> row) -> real_t {
  int nnidx = f(row);
  return sqrt(squareDistBetweenRows(data.GetBlock()[nnidx], row));
}

void run_benchmark(string prefix, function<int(Row<int>)> f,
                   const RowBlockContainer<int> &data, int start) {
  typedef chrono::steady_clock Time;
  typedef chrono::duration<float> fsec;
  real_t total_dist = 0.0;
  auto start_time = Time::now();
  for (int i = start; i < data.Size(); ++i)
    total_dist += query_dist(f, data, data.GetBlock()[i]);
  auto end_time = Time::now();
  real_t mean_dist = total_dist / (data.Size() - start);
  fsec elapsed = end_time - start_time;
  real_t qpersec = real_t(data.Size() - start) / elapsed.count();
  cout << prefix << "mean distance = " << mean_dist << " (" << qpersec
       << " queries per second)." << endl;
}

int main(int argc, char *args[]) {
  int seed = argc > 1 ? atoi(args[1]) : 0;
  int num_to_search = argc > 2 ? atoi(args[2]) : 5000;

  mt19937_64 rng(seed);
  auto data = libsvmread("./mnist.txt");
  vector<int> idxs(num_to_search);
  for (size_t i = 0; i < idxs.size(); ++i)
    idxs[i] = i;
  // Benchmark RPTrees for various values of n0
  for (int n0 = 10; n0 <= num_to_search; n0 *= 2) {
    RandomPartitionTree<int> rpt(rng, 784, n0, data, idxs);
    auto rpt_lookup = bind(&RandomPartitionTree<int>::find_nn, &rpt, _1);
    ostringstream name;
    name << "rpt(n0 = " << n0 << "): ";
    run_benchmark(name.str(), rpt_lookup, data, idxs.size());
  }

  // Benchmark random indexing (mostly for comparing the average distance
  // between points)
  auto random_idx = [&rng, &idxs](Row<int> row) {
    return uniform_int_distribution<int>(0, idxs.size())(rng);
  };
  run_benchmark("random neighbor: ", random_idx, data, idxs.size());
  // Benchmark exhaustive search
  RandomPartitionTree<int> rpt_exhaustive(rng, 784, idxs.size(), data, idxs);
  auto exhaustive_nn =
      bind(&RandomPartitionTree<int>::find_nn, &rpt_exhaustive, _1);
  run_benchmark("exact nn search: ", exhaustive_nn, data, idxs.size());
}
