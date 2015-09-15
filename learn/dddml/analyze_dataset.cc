#include <sstream>
#include <random>

#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"
#include "config_tools.h"


using namespace std;
using namespace dmlc;
using namespace dmlc::data;
using namespace ps;
using namespace dddml;

typedef unsigned long long FeaID;

namespace ps {
  App *App::Create(int argc, char *argv[]) {
    return NULL;
  }
}                               // namespace ps 
/*
 * output function 
 * 
 */ 
void output(std::unordered_map < FeaID, size_t > &counts)
{
  std::multimap < size_t, FeaID > mm;
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    mm.insert(make_pair(it->second, it->first));
  }
  /*
     for (auto it = mm.rbegin(); it != mm.rend(); ++it)
     { 
     cout << it->second << " " << it->first << endl;
     }
   */
}

void output(std::unordered_map < FeaID, size_t > &counts, int maxc,
            std::string featureFile)
{
  std::multimap < size_t, FeaID > mm;
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    mm.insert(make_pair(it->second, it->first));
  }
  int i = 0;
  std::vector < FeaID > topFeatures;
  topFeatures.reserve(maxc);
  for (auto it = mm.rbegin(); it != mm.rend(); ++it) {
    if (i < maxc) {
      topFeatures.push_back(it->second);
      ++i;
    }
    //cout << it->second << " " << it->first << endl;
  }
  dmlc::Stream * file = dmlc::Stream::Create(featureFile.c_str(), "w");
  file->Write(topFeatures);
  delete file;
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "Give conf file as argument. Aborting.. " << std::endl;
    return 0;
  }
  // Read configuration file
  SmartDDDMLConfig cfg = dddml::load_config(argv[1]);
  // Make a random number generator with the given seend
  std::mt19937_64 rng(cfg.safe_analysis_seed());
  // Make a map from features to counts
  std::unordered_map < FeaID, size_t > counts;
  // We will choose parts in the file according to this distribution
  std::uniform_int_distribution < int >part_dist(0,
                                                 cfg.data_parts_per_file() - 1);

  dmlc::data::RowBlockContainer < FeaID > sample;
  dmlc::data::RowBlockContainer < FeaID > *sample_compressed =
    new dmlc::data::RowBlockContainer < FeaID > ();
  unsigned long long count_features = 0;
  for (int fi = 0; fi < cfg.data_num_files(); ++fi) {
    for (int part = 0; part < cfg.analysis_num_parts(); ++part) {
      int partID = part_dist(rng);

      string filename = cfg.data_path(fi);
      std::cout << "Data directory: " << filename << std::endl;
      MinibatchIter < FeaID > reader(filename.c_str(), partID,
                                     cfg.data_parts_per_file(),
                                     cfg.data_format().c_str(),
                                     cfg.analysis_minibatch_size());
      reader.BeforeFirst();
      while (reader.Next()) {
        RowBlock < FeaID > mb = reader.Value();
        dmlc::Localizer < FeaID > lc;
        std::vector < FeaID > uidx;
        std::vector < size_t > freq;
        lc.CountUniqIndex < size_t > (mb, &uidx, &freq);
        for (size_t i = 0; i < uidx.size(); ++i) {
          if (counts.count(uidx[i])) {  // element exists in map
            counts[uidx[i]] += freq[i];
          } else {              // element does not exist
            ++count_features;
            counts[uidx[i]] = freq[i];
          }
        }
      }
    }
  }
  std::cout << "Found " << count_features << " unique features" << std::endl;
  output(counts, cfg.analysis_num_features(), cfg.dispatch_features_path());
  std::cout << "Analysis complete! Output " << cfg.analysis_num_features() <<
    " features." << std::endl;
}
