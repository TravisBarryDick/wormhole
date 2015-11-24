#include "dmlc/data.h"
#include "sample_helper.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include "../base/localizer.h"
#include "../base/minibatch_iter.h"
#include "ps.h"
#include "BufferedWriter.h"
#include <string>

using FeaID = unsigned long long;
using namespace dmlc;
using namespace dmlc::data;

/*
*	Sub-sample data
*
*/


template <typename I>
inline void output(dmlc::ostream* out, RowBlockContainer<I> &rbc)
{
    //std::cout << "Calling output" << std::endl;
    for (int i = 0; i < rbc.Size(); ++i)
    {
      *out << rbc.label[i] ; //label
      for (int j = rbc.offset[i]; j < rbc.offset[i+1]; ++j)
      {
        *out << ' ' << rbc.index[j] ;
        //std::cout << rbc.value << std::endl;
        //std::cout << rbc.index.size() << ' ' << rbc.value.size() << std::endl;
        if (rbc.value.size() != 0) *out  << ':' << rbc.value[j];
      }
      *out << '\n';
    }
}

void translate(std::string inputfile, std::string outputfile)
{
    dmlc::Stream *in = dmlc::Stream::Create(inputfile.c_str(), "r");
    auto out = new dmlc::ostream(dmlc::Stream::Create(outputfile.c_str(), "w"));
    RowBlockContainer<unsigned long long> rbc;
    int len = 0;
    std::cout << len << std::endl;
    while ((rbc.Load(in)))
    {
      std::cout<<"*";
      len += rbc.Size();
      output(out, rbc);
      rbc.Clear();
    }
    std::cout << len << std::endl;
    delete in;
    delete out;
}



void myfunc(
	const char* data_directory,
	const char* outputFile)
{
using real_t = dmlc::real_t;
using namespace dddml;

	/* Step 1: Figure out number of files */
	int nFiles = 1,
		nPartPerFile = 1,
		nPartToRead = 1,
		mb_size = 100000,
		partID;

  char data_format[] = "libsvm";

	int nread = 0, naccept = 0;

	dmlc::data::RowBlockContainer<FeaID> sample;
  BufferedWriter<unsigned long long> writer(outputFile);
  writer.Flush();
	for (int fi = 0; fi < nFiles; ++fi)
	{
		for (int part = 0; part < nPartPerFile; ++part)
		{
			partID = part ;
			//TODO: verify filename
			char filename[200];
			std::sprintf(filename, "%s", data_directory);
			//std::sprintf(filename, "%s/%d", data_directory, fi);
			//std::cout << "fn: " << filename << std::endl;
			//std::sprintf(filename, "%s", data_directory);
			MinibatchIter<FeaID> reader(
				filename, partID, nPartPerFile,
				data_format, mb_size);
			reader.BeforeFirst();
			while (reader.Next()) {
				auto mb = reader.Value(); //row block
				for (size_t i = 0; i < mb.size; ++i)
				{
					++nread;
					sample.Push(mb[i]);
          std::cout << mb.value << std::endl;
          writer.Write(mb[i]);
				}
			}
		}
	}
	//dmlc::Stream *output = dmlc::Stream::Create(outputFile, "w");
	for (int i = 0; i < 100; ++i)
	{
		std::cout << sample.label[i] << " ";
	}
	std::cout << std::endl;

  std::cout << "Done saving " << nread << std::endl;
  


}




namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps

int main(int argc, char *argv[]) {
  //args: i experiment #
  //args: n number of clusters
  //args: train/test
  //char data_dir[] =  "temp1"; //"/home/vpillutl/wormhole/wormhole/learn/data/ctra/train";
  char *dataset = argv[1];
  int i = std::atoi(argv[2]);
  int n = std::atoi(argv[3]);
  char *train = argv[4];
  for (int j = 0; j < n; ++j)
  {
    char input[100]; char output[100];
    sprintf(input, "/home/vpillutl/wormhole/wormhole/learn/data/%s/experiment%d/dispatched/%s/%d/0", dataset, i, train, j);
    sprintf(output, "/home/vpillutl/wormhole/wormhole/learn/data/%s/experiment%d/dispatched/%s/%d_train.txt", dataset, i, train, j);
    //char input[] = "/home/vpillutl/wormhole/wormhole/learn/data/ctra/experiment2/dispatched/train/3/0";
    //char output[] = "temp";
    translate(std::string(input), std::string(output));
  }
}

