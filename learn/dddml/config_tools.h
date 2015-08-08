#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <random>

#include "base/arg_parser.h"
#include "dddml_config.pb.h"

namespace dddml {

namespace {
std::string strjoin(const std::string& a, const std::string& b) {
  std::stringstream buffer;
  buffer << a << b;
  return buffer.str();
}

int make_random_seed() {
  std::random_device rd;
  return rd();
}
};

class SmartDDDMLConfig : public dddmlConfig {
 public:
  std::string dispatch_features_path() const {
    return strjoin(experiment_directory(), dispatch_features_file());
  }

  std::string dispatch_sample_path() const {
    return strjoin(experiment_directory(), dispatch_sample_file());
  }

  std::string dispatch_assignments_path() const {
    return strjoin(experiment_directory(), dispatch_assignments_file());
  }

  std::string dispatch_rpt_path() const {
    return strjoin(experiment_directory(), dispatch_rpt_file());
  }

  std::string dispatched_path() const {
    return strjoin(experiment_directory(), dispatched_directory());
  }

  std::string dispatched_num_clusters_path() const {
    return strjoin(experiment_directory(), dispatched_num_clusters_file());
  }

  int safe_analysis_seed() {
    if (!has_analysis_seed()) set_analysis_seed(make_random_seed());
    return analysis_seed();
  }

  int safe_datasplit_seed() {
    if (!has_datasplit_seed()) set_datasplit_seed(make_random_seed());
    return datasplit_seed();
  }

  int safe_clustering_seed() {
    if (!has_clustering_seed()) set_clustering_seed(make_random_seed());
    return clustering_seed();
  }

  std::string get_data_filename(int part) {
    std::stringstream buffer;
    buffer << data_directory() << part;
    return buffer.str();
  }
};

SmartDDDMLConfig load_config(const std::string& filename) {
  dmlc::ArgParser parser;
  parser.ReadFile(filename.c_str());
  SmartDDDMLConfig cfg;
  parser.ParseToProto(&cfg);
  return cfg;
}
};
