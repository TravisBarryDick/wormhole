#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <random>

#include "base/arg_parser.h"
#include "dddml_config.pb.h"

namespace dddml {
class ConfigWrapper {
 public:
  ConfigWrapper(std::string filename) {
    dmlc::ArgParser parser;
    parser.ReadFile(filename.c_str());
    parser.ParseToProto(&config);
    std::random_device rd;

    data_directory = config.data_directory();
    data_num_files = config.data_num_files();
    data_parts_per_file = config.data_parts_per_file();
    data_num_instances = config.data_num_instances();
    data_dimension = config.data_dimension();
    data_format = config.data_format();

    // paths
    std::string ed = config.experiment_directory();
    dispatch_features_file = join(ed, config.dispatch_features_file());
    dispatch_sample_file = join(ed, config.dispatch_sample_file());
    dispatch_assignments_file = join(ed, config.dispatch_assignments_file());
    dispatch_rpt_file = join(ed, config.dispatch_rpt_file());

    dispatched_directory = join(ed, config.dispatched_directory());
    dispatched_num_clusters_file =
        join(ed, config.dispatched_num_clusters_file());

    // analysis_minibatch_size
    if (config.has_analysis_minibatch_size()) {
      analysis_minibatch_size = config.analysis_minibatch_size();
    } else if (config.has_default_minibatch_size()) {
      analysis_minibatch_size = config.default_minibatch_size();
    } else {
      std::cerr << "Need to specify analysis_minibatch_size" << std::endl;
    }

    // analysis_num_parts
    if (config.has_analysis_num_parts()) {
      analysis_num_parts = config.analysis_num_parts();
    } else if (config.has_default_num_parts()) {
      analysis_num_parts = config.default_num_parts();
    } else {
      std::cerr << "Need to specify analysis_num_parts" << std::endl;
    }

    // analysis_num_features
    if (!config.has_analysis_num_features()) {
      std::cerr << "Need to specify analysis_num_features" << std::endl;
    }
    analysis_num_features = config.analysis_num_features();

    // analysis_seed
    if (config.has_analysis_seed()) {
      analysis_seed = config.analysis_seed();
    } else if (config.has_default_seed()) {
      analysis_seed = config.default_seed();
    } else {
      analysis_seed = rd();
    }

    // datasplit_minibatch_size
    if (config.has_datasplit_minibatch_size()) {
      datasplit_minibatch_size = config.datasplit_minibatch_size();
    } else if (config.has_default_minibatch_size()) {
      datasplit_minibatch_size = config.default_minibatch_size();
    } else {
      std::cerr << "Need to specify datasplit_minibatch_size" << std::endl;
    }

    // datasplit_sample_size
    if (!config.has_datasplit_sample_size()) {
      std::cerr << "Need to specify datasplit_sample_size" << std::endl;
    }
    datasplit_sample_size = config.datasplit_sample_size();

    // datasplit_num_parts
    if (config.has_datasplit_num_parts()) {
      datasplit_num_parts = config.datasplit_num_parts();
    } else if (config.has_default_num_parts()) {
      datasplit_num_parts = config.default_num_parts();
    } else {
      std::cerr << "Need to specify datasplit_num_parts" << std::endl;
    }

    // datasplit_seed
    if (config.has_datasplit_seed()) {
      datasplit_seed = config.datasplit_seed();
    } else if (config.has_default_seed()) {
      datasplit_seed = config.default_seed();
    } else {
      datasplit_seed = rd();
    }

    // clustering_minibatch_size
    if (config.has_clustering_minibatch_size()) {
      clustering_minibatch_size = config.clustering_minibatch_size();
    } else if (config.has_default_minibatch_size()) {
      clustering_minibatch_size = config.default_minibatch_size();
    } else {
      std::cerr << "Need to specify clustering_minibatch_size" << std::endl;
    }

    // clustering_num_clusters
    if (!config.has_clustering_num_clusters()) {
      std::cerr << "Need to specify clustering_num_clusters" << std::endl;
    }
    clustering_num_clusters = config.clustering_num_clusters();

    // clustering_replication
    if (!config.has_clustering_replication()) {
      std::cerr << "Need to specify clustering_replication" << std::endl;
    }
    clustering_replication = config.clustering_replication();

    // clustering_lower_capacity
    if (!config.has_clustering_lower_capacity()) {
      std::cerr << "Need to specify clustering_lower_capacity" << std::endl;
    }
    clustering_lower_capacity = config.clustering_lower_capacity();

    // clustering_upper_capacity
    if (!config.has_clustering_upper_capacity()) {
      std::cerr << "Need to specify clustering_upper_capacity" << std::endl;
    }
    clustering_upper_capacity = config.clustering_upper_capacity();

    // clustering_seed
    if (config.has_clustering_seed()) {
      clustering_seed = config.clustering_seed();
    } else if (config.has_default_seed()) {
      clustering_seed = config.default_seed();
    } else {
      clustering_seed = rd();
    }

    // dispatch_minibatch_size
    if (config.has_dispatch_minibatch_size()) {
      dispatch_minibatch_size = config.dispatch_minibatch_size();
    } else if (config.has_default_minibatch_size()) {
      dispatch_minibatch_size = config.default_minibatch_size();
    } else {
      std::cerr << "Need to specify dispatch_minibatch_size" << std::endl;
    }

    // dispatch_rpt_n0
    if (!config.has_dispatch_rpt_n0()) {
      std::cerr << "Need to specify dispatch_rpt_n0" << std::endl;
    }
    dispatch_rpt_n0 = config.dispatch_rpt_n0();
  }

  /** returns a string for the `filenum`th file for the data */
  std::string get_data_filename(int filenum) {
    std::stringstream buffer;
    buffer << data_directory << filenum;
    return buffer.str();
  }

  std::string data_directory;
  int data_num_files;
  int data_parts_per_file;
  int data_num_instances;
  int data_dimension;
  std::string data_format;

  std::string dispatch_features_file;
  std::string dispatch_sample_file;
  std::string dispatch_assignments_file;
  std::string dispatch_rpt_file;

  std::string dispatched_directory;
  std::string dispatched_num_clusters_file;

  int analysis_minibatch_size;
  int analysis_num_parts;
  int analysis_num_features;
  int analysis_seed;

  int datasplit_minibatch_size;
  int datasplit_sample_size;
  int datasplit_num_parts;
  int datasplit_seed;

  int clustering_minibatch_size;
  int clustering_num_clusters;
  int clustering_replication;
  float clustering_lower_capacity;
  float clustering_upper_capacity;
  int clustering_seed;

  int dispatch_minibatch_size;
  int dispatch_rpt_n0;

 private:
  std::string join(std::string a, std::string b) {
    std::stringstream buffer;
    buffer << a << b;
    return buffer.str();
  }

  dddml::dddmlConfig config;
};
}
