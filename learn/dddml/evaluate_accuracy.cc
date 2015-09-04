#include "config_tools.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <dmlc/io.h>
#include "data/row_block.h"
using FeaID = unsigned long long;


int main(int argc, char *argv[])
{
  if (argc < 4) {
    std::cerr << "The Correct usage is: " << argv[0] <<
      " <config_file>  <num_machines> <num_parts>" << std::endl;
    return 0;
  }
  dddml::SmartDDDMLConfig cfg = dddml::load_config(argv[1]);
  std::istringstream iss(argv[2]);
  int num_parts, num_machines;
  bool temp = (iss >> num_machines);
  if (!temp) {
    std::cerr << "The Correct usage is: " << argv[0] <<
      " <config_file>  <num_machines> <num_parts>" << std::endl;
    return 0;
  }
  std::istringstream iss1(argv[3]);
  temp = (iss1 >> num_parts);
  if (!temp) {
    std::cerr << "The Correct usage is: " << argv[0] <<
      " <config_file>  <num_machines> <num_parts>" << std::endl;
    return 0;
  }
  //eval accuracy
  //filename: $(machine_id)/$(part_id)_part-0
  unsigned long long total = 0, correct = 0;

  for (int i = 0; i < num_machines; ++i) {
    for (int j = 0; j < num_parts; ++j) {
      //step 1: filenames
      std::ostringstream oss, oss1;
      oss << cfg.predictions_path() << i << "/" << j << "_part-0";
      auto pred_file = oss.str();
      std::cerr << "**pred_path: " << pred_file << std::endl;
      oss1 << cfg.dispatched_path(i, true) << j;        // is_test=true gives testing path
      auto val_file = oss1.str();
      std::cerr << "val_file: " << val_file << std::endl;
      //step 2: read true labels
      dmlc::Stream * out = dmlc::Stream::Create(val_file.c_str(), "r");
      dmlc::data::RowBlockContainer < FeaID > rbc;
      rbc.Load(out);
      delete out;
      //step 3: read predictions and evaluate accuracy
      int index = 0;
      std::ifstream predfile(pred_file);
      for (double a; predfile >> a;) {
        //labels: +1 or -1
        //if (cfg.prob_predict()) a -= 0.5; // convert prob to > 0 and < 0.
        ++total;
        if (rbc.label[index] * a > 0) {
          ++correct;
        }
        ++index;
      }
    }
  }
  double accuracy =
    static_cast < double >(correct) / static_cast < double >(total);
  std::cerr << "Final Accuracy: " << accuracy << std::endl;
  std::cout << accuracy << std::endl;
}
