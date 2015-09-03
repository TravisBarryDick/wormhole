# dddml

To do before running:
- data_directory of the conf file has the data
	- $(data_directory)/train contains training files numbered 0, 1,...
	- $(data_directory)/test contains the testing files
	- number of training files is set in field "datasplit_num_parts" of the conf file
- experiment_directory exists: it will contain all the temporaries created
- $(experiment_directory)/model, $(experiment_directory)/dispatched, $(e)/dispatched/train, $(e)/dispatched/test, $(e)/dispatched/train/$(i), $(e)/dispatch/test/$(i) should be created too 

