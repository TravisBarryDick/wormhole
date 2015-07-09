#include <sstream>

#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "base/arg_parser.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"

#include "dispatcher_config.pb.h"

using namespace std;
using namespace dmlc;
using namespace dmlc::data;
using namespace ps;
using namespace dddml;

typedef unsigned FeaID;

int CreateServerNode(int argc, char* argv[]) { return 0; }

int WorkerNodeMain(int argc, char* argv[]) {
  ArgParser parser;
  if (argc > 1 && strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
  parser.ReadArgs(argc - 2, argv + 2);
  DispatcherConfig conf;
  parser.ParseToProto(&conf);

  // Steps:
  // 1. Read the dispatch data, cluster assignment, and tree
  // 2. Use MyRank() and RankSize() to decide what parts of the file to read
  // 3. Iterate through the part belonging to this worker and dispatch.

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

  for (int file_num = 0; file_num < conf.n_files(); ++file_num) {
    stringstream filename;
    filename << conf.input_data() << file_num;
    cout << "Worker " << MyRank() << " opening parts " << first_part << "-"
         << (first_part + num_parts - 1) << " of " << filename.str() << endl;
    for (size_t part = first_part; part < first_part + num_parts; ++part) {
      MinibatchIter<FeaID> reader(filename.str().c_str(), part,
                                  static_cast<size_t>(conf.n_parts_per_file()),
                                  conf.input_format().c_str(),
                                  static_cast<size_t>(conf.mb_size()));
      reader.BeforeFirst();
      while (reader.Next()) {
        auto mb = reader.Value();  // mb is a RowBlock        
      }
    }
  }
  cout << "Worker " << MyRank() << " done!" << endl;
}
