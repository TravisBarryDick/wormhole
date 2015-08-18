#include "config_tools.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include "ps.h"

void create_training_conf(std::string config_filename, std::string train_data, std::string model_file, int machine_id){
  //write these to file
  //TODO:check double quotes of output:
  std::ofstream out1(config_filename, std::ofstream::out);
  out1 << "train_data = \"" << train_data << "\"\n"; 
  out1 << "model_out = \"" << model_file << machine_id << "\"\n";
  out1 << "data_format = \"RowBlockContainer\"\n";
  out1.close();
}

void create_testing_conf(std::string config_filename, std::string val_data, std::string pred_file, std::string model_file, int machine_id){
  //write these to file
  //TODO:check double quotes of output:
  std::ofstream out1(config_filename, std::ofstream::out);
  out1 << "val_data = \"" << val_data << "\"\n"; 
  out1 << "model_in = \"" << model_file << machine_id << "\"\n";
  out1 << "data_format = \"RowBlockContainer\"\n";
  out1 << "pred_out = \"" << pred_file << "\"\n";
  out1.close();
}


namespace ps {
App* App::Create(int argc, char *argv[]) {
  return NULL;
}
}  // namespace ps



int main(int argc, char *argv[]) {
  if (argc < 4){
  	std::cout << "The Correct usage is: " << argv[0] << " <config_file> <train_or_test> <machine_id> [num_parts{for test only}]" << std::endl;
  	return 0;
  }
  dddml::SmartDDDMLConfig cfg = dddml::load_config(argv[1]);
  std::string train_or_test (argv[2]);
  std::istringstream iss( argv[3] );
  int machine_id;
  bool temp = (iss >> machine_id);
  if (!temp) 
  {
  	std::cout << "The Correct usage is: " << argv[0] << "<config_file> <train_or_test> <machine_id> [num_parts{for test only}]" << std::endl;
  	return 0;
  }

  std::string conf_filename ("temp_LEARN.conf");
  std::string script_command ("tracker/dmlc_local.py"); //TODO
  std::stringstream sss;
  sss << script_command << " -n 1 -s 1 bin.linear.dmlc " << conf_filename;
  std::string command = sss.str();
  //std::string command ("tracker/dmlc_local.py -n 1 -s 1 bin/linear.dmlc <conf_file>"); //TODO
  
  if (train_or_test.compare("train") == 0)
  {
  	std::cout << "training file" << std::endl;
  	//create conf file
  	create_training_conf(conf_filename, cfg.dispatched_path(machine_id, false), cfg.model_path(), machine_id);
  	//run command
  	std::system(command.c_str()); 
  	//std::cout << (command) << std::endl; 
  	
  }
  else
  {
  	std::cout << "testing file" << std::endl;
  	if (argc < 5) 
    {
    	std::cout << "The Correct usage is: " << argv[0] << "<config_file> <train_or_test> <machine_id> [num_parts{for test only}]" << std::endl;
    	return 0;
    }
  	std::istringstream iss1( argv[4] );
    int num_parts;
    temp = (iss1 >> num_parts);
    if (!temp) 
    {
    	std::cout << "The Correct usage is: " << argv[0] << "<config_file> <train_or_test> <machine_id> [num_parts{for test only}]" << std::endl;
    	return 0;
    }
    
    for (int i = 0; i < num_parts; ++i)
    {
      std::ostringstream oss, oss1;
      oss << cfg.predictions_path() << i;
      auto pred_directory = oss.str();
      oss1 << cfg.dispatched_path(machine_id, true) << i;
      auto val_file = oss1.str();
    	//create conf file
    	create_testing_conf(conf_filename, val_file, pred_directory, cfg.model_path(), machine_id);
    	//run command
  		std::system(command.c_str()); 
    	//std::cout << (command) << std::endl; 
  	}
  }
  
}


