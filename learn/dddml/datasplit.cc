
#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"
#include "base/arg_parser.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <random>
#include <vector>
#include "sample_helper.h"

#include "dddml_config.pb.h"

namespace dddml{
using FeaID = unsigned long long;
using namespace dmlc;
using namespace dmlc::data;

/*
*	Sub-sample data
*
*/


void ReadFile(const char* featureFile, std::vector<FeaID> *features)
{
	dmlc::Stream *file = dmlc::Stream::Create(featureFile, "r", true);
	if (file == NULL)
	{
		std::cerr << "Feature File: " << featureFile << ", doesn't exist\n";
		//exit(-1);
	}
	else
	{
		file->Read(features);
		std::cout << "DEBUG: "; for (auto pvix : *features) std::cout << pvix << " "; std::cout << std::endl;
	}
	delete file;
}

std::vector<FeaID> *Intersect(std::vector<FeaID> &v1, std::vector<FeaID> &v2)
{
	// features do not exist
	if (v1.size() == 0) return &v2;
	else if (v2.size() == 0) return &v1;
	//else:
	std::vector<FeaID> *output = new std::vector<FeaID>();
	for (unsigned i=0,j=0;((i < v1.size()) && (j < v2.size())); )
	{
		if (v1[i] == v2[j])
		{
			output->push_back(v1[i]);
			++i; ++j;
		}
		else if(v1[i] > v2[j])
		{
			++j;
		}
		else if (v1[i] < v2[j])
		{
			++i;
		}
	}
	return output;
}

void subsample(
	const char* featureFile, //name of feature file
	const char* data_directory,
	const char* outputFile,
	const char *data_format,
	unsigned int subsample_size,
	unsigned int total_size,
	std::mt19937_64 &rng,
	int nFiles,
	int nPartPerFile,
	int nPartToRead,
	int mb_size
)
{
using real_t = dmlc::real_t;

	/* Step 1: Figure out number of files */
	int /*nFiles = 1,
		nPartPerFile = 100,
		nPartToRead = 10,
		mb_size = 1000,*/
		partID;

	std::uniform_real_distribution<> dist(0, 1);
	std::uniform_int_distribution<int> dis(0, nPartPerFile - 1);

	real_t probability_of_selecting_one_row = (static_cast<real_t> (subsample_size)) / total_size * nPartPerFile / nPartToRead; //TODO: Check.
	std::cerr << "Probability of selecting a row: " << probability_of_selecting_one_row << std::endl;
	probability_of_selecting_one_row = (probability_of_selecting_one_row > 1.0) ? 1.0 : probability_of_selecting_one_row;

	/* Step 2: Read some of the blocks at random, and sub-sample */

	int nread = 0, naccept = 0;

	dmlc::data::RowBlockContainer<FeaID> sample;
	dmlc::data::RowBlockContainer<FeaID> *sample_compressed = new dmlc::data::RowBlockContainer<FeaID>();

	for (int fi = 0; fi < nFiles; ++fi)
	{
		for (int part = 0; part < nPartToRead; ++part)
		{
			std::cout << "DEBUGGING: new part" << std::endl;
			partID = dis(rng);
			//TODO: verify filename
			char filename[200];
			//std::sprintf(filename, "%s/%d", data_directory, fi);
			std::sprintf(filename, "%s", data_directory);
			MinibatchIter<FeaID> reader(
				filename, partID, nPartPerFile,
				data_format, mb_size);
			std::cout << "DEBUGGING: starting read\n";
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				std::cout << "DEBUGGING: read a row block " << std::endl;
				for (size_t i = 0; i < mb.size; ++i)
				{
					//decide whether to add row mb[i] to the sample or not
					++nread;
					if (dist(rng) < probability_of_selecting_one_row)
					{
						++naccept;
						sample.Push(mb[i]);
					}

				}
			}
		}
	}

	/* Step 3: Localize */
	dmlc::RowBlock<FeaID> sample1 = sample.GetBlock();
	/* 3.1: read feature file */
	//int SomeDefaultStartingValue = 10000; //TODO
	std::vector<FeaID> features;
	//features.reserve(SomeDefaultStartingValue);
	ReadFile(featureFile, &features);
	/* 3.2: Get set of features to keep using localizer */
	dmlc::Localizer <FeaID> lc;
	std::vector<FeaID> *uidx = new std::vector<FeaID>();
	lc.CountUniqIndex<FeaID>(sample1, /*4,*/ uidx, NULL); //include nthreads = 4 for older version

	std::cout << "DEBUGGING: features.size() = " << features.size()
						<< " uidx->size() = " << uidx->size() << std::endl;

	/* 3.3: intersect uidx with features */
	std::vector<FeaID> *idx_dict = Intersect(features, *uidx);

	std::cout << "DEBUGGING: ";
	for (auto i : *idx_dict) std::cout << i << " ";
	std::cout << std::endl;

	/* 3.4: localize */
	lc.RemapIndex<FeaID>(sample1, *idx_dict, sample_compressed);


	 std::cout << "DEBUG: sample_compressed->Size() = " << sample_compressed->Size() << std::endl;
	// for (size_t i = 0; i < sample_compressed->Size(); ++i) {
	// 	std::cout << sample_compressed->GetBlock()[i].length << " ";
	// }
	// std::cout << std::endl;

	/* Step 4: Write to file */

	//First write idx_dict. Then write compressed sample.

	dmlc::Stream *output = dmlc::Stream::Create(outputFile, "w");
	output->Write(*idx_dict);
	sample_compressed->Save(output);
	std::cout << "Done sampling" << std::endl;
	delete output;
}

} //namespace dddml



namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps


int main(int argc, char *argv[])
{
	using namespace dddml;
	ArgParser parser;
	if (argc > 1 && strcmp(argv[1], "none")) parser.ReadFile(argv[1]);
	parser.ReadArgs(argc - 2, argv + 2);
	dddmlConfig conf;
	parser.ParseToProto(&conf);

	#if 0
	char featureFile[] = "features.txt";
	char data_directory[] = "../data/mnist";
	char outputFile[] = "../data/mnist.out";
	char data_format[] = "libsvm";
	int subsample_size = 10000;
	int total_size = 60000;
	#endif

	const char *featureFile = conf.feature_filename().c_str();
	const char *data_directory = conf.data_directory().c_str();
	const char *outputFile = conf.sample_filename().c_str();
	const char *data_format = conf.data_format().c_str();
	int subsample_size = conf.sample_size();
	int total_size = conf.total_size_of_dataset();
	//std::cout << subsample_size << " " << total_size << " " << featureFile << std::endl;

	//int nfile = conf.n_files(),

	std::random_device rd;
	int seed = conf.seed();
        std::mt19937_64 rng (seed);

	subsample(featureFile, data_directory, outputFile,data_format, subsample_size, total_size, rng,
			conf.n_files(), conf.n_parts_per_file(), conf.n_parts_to_read(), conf.analysis_minibatch_size());


	return 0;
}
