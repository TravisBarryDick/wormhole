include ../../make/ps_app.mk

CFLAGS  = -O3 -ggdb -std=c++11 -I./ -I../ $(INCLUDE) $(DMLC_CFLAGS) $(PS_CFLAGS) $(EXTRA_CFLAGS) #-Wall
LDFLAGS = $(CORE_PATH)/libdmlc.a $(DMLC_LDFLAGS) $(PS_LDFLAGS) $(EXTRA_LDFLAGS)

all: base core bin ps # tool # test

clean:
	rm -rf build/*


bin: build/datasplit build/analyze_criteo

#build/%.o: %.cc
#	@mkdir -p $(@D)
#	$(CXX) $(CFLAGS) -MM -MT build/$*.o $< >build/$*.d
#	$(CXX) $(CFLAGS) -c $< -o $@

build/%.o: %.cc
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -MM -MT build/$*.o $< >build/$*.d
	$(CXX) $(CFLAGS) -c $< -o $@


%.pb.cc %.pb.h : %.proto
	${DEPS_PATH}/bin/protoc --cpp_out=. --proto_path=. $<

build/datasplit: build/dddml_config.pb.o build/datasplit.o  ../base/base.a   $(PS_PATH)/build/libps.a
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/trial: build/trial.o  ../base/base.a   $(PS_PATH)/build/libps.a
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/dispatcher: build/dispatcher.o build/dddml_config.pb.o ../base/base.a RPTS.h $(PS_PATH)/build/libps.a
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

build/analyze_criteo: build/dddml_config.pb.o build/analyze_criteo.o ../base/base.a RPTS.h $(PS_PATH)/build/libps.a
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@


build/kmeans: build/dddml_config.pb.o build/kmeans.o  ../base/base.a   $(PS_PATH)/build/libps.a
	$(CXX) $(CFLAGS) $(filter %.o %.a, $^) $(LDFLAGS) -o $@

-include build/*.d
