include ../../make/ps_app.mk

CFLAGS  = -O3 -ggdb -std=c++11 -I./ -I../ $(INCLUDE) $(DMLC_CFLAGS) $(PS_CFLAGS) $(EXTRA_CFLAGS) #-Wall
LDFLAGS = $(CORE_PATH)/libdmlc.a $(DMLC_LDFLAGS) $(PS_LDFLAGS) $(EXTRA_LDFLAGS)

all: bin # tool # test

clean:
	rm -rf build/*

bin: build/datasplit build/analyze_criteo build/kmeans build/dispatcher


build/%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -MM -MT build/$*.o $< >build/$*.d
	$(CXX) $(CFLAGS) -c $< -o $@

%.pb.cc %.pb.h : %.proto
	${DEPS_PATH}/bin/protoc --cpp_out=. --proto_path=. $<

build/datasplit: build/dddml_config.pb.o build/datasplit.o  $(DMLC_SLIB)
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/trial: build/dddml_config.pb.o build/trial.o $(DMLC_SLIB) 
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/dispatcher: build/dddml_config.pb.o build/dispatcher.o $(DMLC_SLIB)
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/analyze_criteo: build/dddml_config.pb.o build/analyze_criteo.o $(DMLC_SLIB)
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@


build/kmeans: build/dddml_config.pb.o build/kmeans.o  $(DMLC_SLIB)
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

-include build/*.d
