#include <sstream>

#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"
#include "kmeans.h"
#include "RPTS.h"
#include "BufferedWriter.h"
#include "config_tools.h"

using namespace std;
using namespace dmlc;
using namespace dmlc::data;
using namespace ps;
using namespace dddml;

typedef unsigned long long FeaID;

int CreateServerNode(int argc, char* argv[]) { return 0; }

int WorkerNodeMain(int argc, char* argv[]) {
  ConfigWrapper cfg(argv[1]);

  // Load the dispatching data and the dispatch tree
  RandomPartitionTree<FeaID> rpt(cfg.dispatch_sample_file.c_str(),
                                 cfg.dispatch_rpt_file.c_str());

  // Load the cluster assignments for the dispatch data
  vector<vector<int>> assignments;
  read_assignments(cfg.dispatch_assignments_file.c_str(), &assignments);

  int parts_per_worker = cfg.data_parts_per_file / RankSize();
  int extra_parts = cfg.data_parts_per_file - RankSize() * parts_per_worker;

  // The number of parts per file that this worker is responsible for
  int num_parts = parts_per_worker;
  if (MyRank() < extra_parts) {
    num_parts += 1;
  }
  // The first part per file that this worker is responsible for
  int first_part = MyRank() * parts_per_worker;
  first_part += min(extra_parts, MyRank());

  // Read the number of clusters
  int k;
  ifstream num_clusters(cfg.dispatched_num_clusters_file, ios::in);
  num_clusters >> k;
  num_clusters.close();

  // Make a BufferedWriter for each cluster
  vector<BufferedWriter<FeaID>> cluster_writers;
  cluster_writers.reserve(k);
  for (size_t i = 0; i < k; ++i) {
    stringstream filename;
    filename << cfg.dispatched_directory << i << "/" << MyRank();
    string filename_str = filename.str();
    cluster_writers.push_back(BufferedWriter<FeaID>(filename_str.c_str()));
  }

  for (int file_num = 0; file_num < cfg.data_num_files; ++file_num) {
    string filename = cfg.get_data_filename(file_num);
    for (size_t part = first_part; part < first_part + num_parts; ++part) {
      MinibatchIter<FeaID> reader(
          filename.c_str(), part, static_cast<size_t>(cfg.data_parts_per_file),
          cfg.data_format.c_str(), cfg.dispatch_minibatch_size);
      reader.BeforeFirst();
      while (reader.Next()) {
        auto mb = reader.Value();
        for (size_t i = 0; i < mb.size; ++i) {
          size_t nn_idx = rpt.find_nn(mb[i]);
          for (auto cluster_ix : assignments[nn_idx]) {
            cluster_writers[cluster_ix].Write(mb[i]);
          }
        }
      }
    }
  }
  for (auto cw : cluster_writers) cw.Flush();
}
