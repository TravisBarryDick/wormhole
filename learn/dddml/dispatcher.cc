#include <sstream>

#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "base/arg_parser.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"
#include "kmeans_helper.h"
#include "RPTS.h"
#include "BufferedWriter.h"

#include "dddml_config.pb.h"

using namespace std;
using namespace dmlc;
using namespace dmlc::data;
using namespace ps;
using namespace dddml;

typedef unsigned long long FeaID;

int CreateServerNode(int argc, char* argv[]) { return 0; }

int WorkerNodeMain(int argc, char* argv[]) {
  ArgParser parser;
  if (argc > 1 && strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
  parser.ReadArgs(argc - 2, argv + 2);
  dddmlConfig conf;
  parser.ParseToProto(&conf);

  // Load the dispatching data and the dispatch tree
  auto rpt = RandomPartitionTree<FeaID>(conf.sample_filename().c_str(),
                                        conf.rpt_filename().c_str());
  // Load the cluster assignments for the dispatch data
  vector<vector<int>> assignments;
  read_assignments(conf.assignments_filename().c_str(), &assignments);

  int parts_per_worker = conf.n_parts_per_file() / RankSize();
  int extra_parts = conf.n_parts_per_file() - RankSize() * parts_per_worker;

  // The number of parts per file that this worker is responsible for
  int num_parts = parts_per_worker;
  if (MyRank() < extra_parts) {
    num_parts += 1;
  }
  // The first part per file that this worker is responsible for
  int first_part = MyRank() * parts_per_worker;
  first_part += min(extra_parts, MyRank());


  cout << "GOT THIS FAR" << endl;

  // Read the number of clusters
  int k;
  ifstream num_clusters(conf.num_clusters_filename().c_str(), ios::in);
  num_clusters >> k;
  num_clusters.close();

  // Make a BufferedWriter for each cluster
  vector<BufferedWriter<FeaID>> cluster_writers;
  cluster_writers.reserve(k);
  for (size_t i = 0; i < k; ++i) {
    stringstream filename;
    filename << conf.final_output_directory() << i << "/" << MyRank();
    string filename_str = filename.str();
    cluster_writers.push_back(BufferedWriter<FeaID>(filename_str.c_str()));
  }

  for (int file_num = 0; file_num < conf.n_files(); ++file_num) {
    stringstream filename;
    filename << conf.data_directory() << file_num;
    auto filename_str = filename.str();
    for (size_t part = first_part; part < first_part + num_parts; ++part) {
      MinibatchIter<FeaID> reader(
          filename_str.c_str(), part,
          static_cast<size_t>(conf.n_parts_per_file()),
          conf.data_format().c_str(),
          static_cast<size_t>(conf.dispatch_minibatch_size()));

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
}
