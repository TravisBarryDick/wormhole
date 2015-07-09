#include <iostream>
#include <cstring>
#include <algorithm>
#include <random>
#include <vector>

#include "sample_helper.h"
#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "base/arg_parser.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"

#include "dispatcher_config.pb.h"

using namespace std;

class DispatcherWorker : public ps::App {
 public:
  DispatcherWorker(dddml::DispatcherConfig& confg) {}

  virtual bool Run() {
    cout << "Started the worker!\n";
    return true;
  }
};

class DispatcherServer : public ps::App {
 public:
  DispatcherServer(dddml::DispatcherConfig& conf) {}

  virtual bool Run() {
    cout << "Started the server!\n";
    return true;
  }
};

class DispatcherScheduler : public ps::App {
 public:
  DispatcherScheduler(dddml::DispatcherConfig& conf) {}

  virtual bool Run() {
    cout << "Started the scheduler!\n";
    cout << "There are " << ps::NumServers() << " servers and "
         << ps::NumWorkers() << " workers.\n";
    return true;
  }
};

namespace ps {
App* App::Create(int argc, char* argv[]) {
  CHECK_GE(argc, 2) << "\nusage: " << argv[0] << " conf_file\n";
  dmlc::ArgParser parser;
  if (strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
  parser.ReadArgs(argc - 2, argv + 2);
  dddml::DispatcherConfig conf;
  parser.ParseToProto(&conf);
  if (IsWorkerNode()) {
    return new DispatcherWorker(conf);
  } else if (IsServerNode()) {
    return new DispatcherServer(conf);
  } else if (IsSchedulerNode()) {
    return new DispatcherScheduler(conf);
  }
};
};

int main(int argc, char* argv[]) { return ps::RunSystem(&argc, &argv); }
