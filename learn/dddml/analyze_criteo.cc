#include <sstream>
#include <random>

#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "base/arg_parser.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"

#include "dddml_config.pb.h"

using namespace std;
using namespace dmlc;
using namespace dmlc::data;
using namespace ps;
using namespace dddml;

typedef unsigned FeaID;

namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps


void output (std::unordered_map<FeaID, size_t> &counts)
{
  std::multimap<size_t, FeaID> mm;
  for (auto it = counts.begin(); it != counts.end(); ++it)
  {
    mm.insert(make_pair(it->second, it->first));
  }
  for (auto it = mm.rbegin(); it != mm.rend(); ++it)
	{
	  cout << it->second << " " << it->first << endl;
	}
}

void output(std::unordered_map<FeaID, size_t> &counts, int maxc, std::string &featureFile)
{
  std::multimap<size_t, FeaID> mm;
  for (auto it = counts.begin(); it != counts.end(); ++it)
  {
    mm.insert(make_pair(it->second, it->first));
  }
  int i = 0;
  std::vector<FeaID> topFeatures; topFeatures.reserve(maxc);
  for (auto it = mm.rbegin(); it != mm.rend(); ++it)
	{
	  if (i <= maxc)
	  {
	    topFeatures.push_back(it->second);
	    ++i;
	  }
	  cout << it->second << " " << it->first << endl;
	}
	dmlc::Stream *file = dmlc::Stream::Create(featureFile.c_str(), "w");
	file->Write(topFeatures);
}

int main(int argc, char *argv[])
{
  ArgParser parser;
  if (argc > 1 && strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
  parser.ReadArgs(argc - 2, argv + 2);
  dddmlConfig conf;
  parser.ParseToProto(&conf);

  std::random_device rd;
	std::mt19937_64 rng (rd());

  // Steps:
  // 1. Read the dispatch data, cluster assignment, and tree
  // 2. Use MyRank() and RankSize() to decide what parts of the file to read
  // 3. Iterate through the part belonging to this worker and dispatch.

  int nFiles = conf.n_files(),
		nPartPerFile = conf.n_parts_per_file(),
		nPartToRead = conf.n_parts_to_read(),
		mb_size = conf.analysis_minibatch_size(),
		n_features = conf.n_features_to_pick(),
		partID;
  auto data_directory = conf.data_directory();
  auto data_format = conf.data_format();
  auto feature_filename = conf.feature_filename();

  std::unordered_map<FeaID, size_t> counts;

	std::uniform_int_distribution<int> dis(0, nPartPerFile - 1);

	real_t probability_of_selecting_one_row = (static_cast<real_t> (nPartToRead) / nPartPerFile);
	std::cerr << "Probability of selecting a row: " << probability_of_selecting_one_row << std::endl;
	probability_of_selecting_one_row = (probability_of_selecting_one_row > 1.0) ? 1.0 : probability_of_selecting_one_row;

	dmlc::data::RowBlockContainer<FeaID> sample;
	dmlc::data::RowBlockContainer<FeaID> *sample_compressed = new dmlc::data::RowBlockContainer<FeaID>();

	for (int fi = 0; fi < nFiles; ++fi)
	{
		for (int part = 0; part < nPartToRead; ++part)
		{
			partID = dis(rng);
			//TODO: verify filename
      stringstream filename;
      filename << data_directory << fi;
      auto filename_str = filename.str();
			MinibatchIter<FeaID> reader(
				filename_str.c_str(), partID, nPartPerFile,
				data_format.c_str(), mb_size);
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				dmlc::Localizer <FeaID> lc;
	      std::vector<FeaID> uidx;
	      std::vector<size_t> freq;
	      lc.CountUniqIndex<size_t>(mb, &uidx, &freq);
	      for (size_t i = 0; i < uidx.size(); ++i)
	      {
	        if (counts.count(uidx[i])) {// element exists in map
	          counts[uidx[i]] += freq[i];
	        }
	        else{ //element does not exist
	          counts[uidx[i]] = freq[i];
	        }
	      }
			}
		}
	}

	/*
	for (auto it: counts)
	{
	  cout << it.first << " " << it.second << endl;
	}
	*/
	output(counts, n_features, feature_filename);

}
