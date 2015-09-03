#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
//#include "ps.h"
#include <limits>


/* TODO:
- Multi-class learning:
- Give training and testing configuration files, but without train_data and val_data specified.
- Pass original training and testing file names

Algo:
produce conf files
foreach (i in range)
	produce new training and testing files
	run train and test to get pred_i
combine all pred_i to get one predictions file
compute accuracy versus original testing file

*/

void create_training_config(std::string config_filename, std::string old_config_file, std::string train_filename){
  //write these to file
  std::ofstream out1(config_filename, std::ofstream::out);
  out1 << "train_data = \"" << train_filename << "\"\n";
  
  std::ifstream in1(old_config_file);
  for( std::string line; getline( in1, line ); )
  {
    out1 << line << "\n";
  }
  in1.close();
  out1.close();
}


void create_testing_config(std::string config_filename, std::string old_config_file, std::string test_filename, std::string pred_filename){
  //write these to file
  std::ofstream out1(config_filename, std::ofstream::out);
  out1 << "val_data = \"" << test_filename << "\"\n";
  out1 << "predict_out = \"" << pred_filename << "\"\n";
  std::ifstream in1(old_config_file);
  for( std::string line; getline( in1, line ); )
  {
    out1 << line << "\n";
  }
  
  in1.close();
  out1.close();
}

void convert_to_binary_class(std::string original_filename, std::string new_filename, int positive_class)
{
  std::ifstream in1(original_filename, std::ifstream::in);
  std::ofstream out1(new_filename, std::ofstream::out);
  for( std::string line; getline( in1, line ); )
  {
    int digit = line.at(0) - '0'; //class is single digit
    line.at(0) = (digit == positive_class) ? '1' : '0'; // for 1 or 0; new class label
    out1 << line << "\n";
  }
  in1.close();
  out1.close();

}

void combine_predictions(std::string pred_file, std::string final_predictions)
{
  int nclass = 10;
  std::ifstream in[nclass];
  for(int i = 0; i < nclass; i++){
    std::stringstream filename;
    filename << pred_file << i ;
    in[i].open(filename.str());
  }
  std::ofstream out(final_predictions, std::ofstream::out);
  
  bool reading = true;
  std::string line;
  while (reading)
  {
  int max_index = 0;
  float max_val = -1 * std::numeric_limits<float>::max();
  for (int i = 0; i < nclass; ++i)
  {
    reading = reading && getline( in[i], line);
    if (!reading) break;
    float pred = std::stof(line);
    if (max_val < pred) // output
    {
      max_val = pred;
      max_index = i;
    }
  }
  out << max_index << "\n";
  }
  out.close();
  for (int i = 0; i < nclass; ++i) in[i].close(); 
}


float eval_accuracy(std::string final_predictions, std::string test_file)
{
  std::ifstream in1(final_predictions, std::ifstream::in);
  std::ifstream in2(test_file, std::ifstream::in);
  
  std::string line1, line2;
  unsigned long long total = 0, correct = 0;
  while(getline(in1, line1) && getline(in2, line2))
  {
    ++total;
    if (line1.at(0) == line2.at(0)) ++correct;
  }
  in1.close(); in2.close();
  return ((float) correct) / total;
}


/*
namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps
*/

int main(int argc, char *argv[]) {
  if (argc < 6){
  	std::cout << "The Correct usage is: " << argv[0] << " <train_config_file> <test_config_file> <train_data> <val_data> <linear_or_difacto>" << std::endl;
  	return 0;
  }
  std::string train_config_file (argv[1]);
  std::string test_config_file (argv[2]);
  std::string train_file (argv[3]);
  std::string test_file (argv[4]);
  std::string linear_or_difacto (argv[5]);
  
  std::string new_train_conf ("temp_train.conf");
  std::string new_test_conf ("temp_test.conf");
  
  std::string new_train_file("train.txt");
  //std::string new_test_file("test.txt");
  
  //TODO: fix this
  std::string script_command ("../../tracker/dmlc_local.py");
  std::stringstream sss, sss1, sss2;
  sss << script_command << " -n 1 -s 1 ../../bin/" << linear_or_difacto << ".dmlc ";
  sss1 << sss.str() << new_train_conf;
  sss2 << sss.str() << new_test_conf;
  std::string train_command = sss1.str(), test_command = sss2.str();
  std::string pred_file("pred_");
  
  //create train config
  create_training_config(new_train_conf, train_config_file, new_train_file);
  for(int i = 0; i < 10; ++i)
  {
    std::cout << "****************\nStarting " << i << "\n*******************" << std::endl;
    std::system("rm model*");
    //create new train and test files
    convert_to_binary_class(train_file, new_train_file, i);
    //convert_to_binary_class(test_file, new_test_file, i);
    
    //create new test config
    std::stringstream ss;
    ss << pred_file << i;
    create_testing_config(new_test_conf, test_config_file, /*new_*/test_file, ss.str());
    
    //train and test
    std::cout << train_command << std::endl;
    std::cout << test_command << std::endl;
    std::system(train_command.c_str());
    std::system(test_command.c_str());
    char cmd[200];
    std::sprintf(cmd, "cat pred_%dt* > pred_%d; rm pred_%dt*", i, i,i);
    std::system(cmd);
    //return 0;
  }
  //combine predictions
  std::string final_predictions = (pred_file + "final");
  combine_predictions(pred_file, final_predictions);
  
  //evaluate accuracy
  float accuracy = eval_accuracy(final_predictions, test_file);  
  std::cout << "accuracy: " << accuracy << std::endl;  
}


