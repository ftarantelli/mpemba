LIBRARY_DIR=$(HOME)/Desktop/C++/itensor
TDVP_DIR=$(LIBRARY_DIR)/tdvp

APP=mpsising

HEADERS=$(TDVP_DIR)/tdvp.h $(PWD)/$(APP).hpp

CCFILES=$(APP).cpp

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

CCFLAGS+=-I$(TDVP_DIR)
CCGFLAGS+=-I$(TDVP_DIR)

#Mappings --------------
OBJECTS=$(patsubst %.cpp,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cpp $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cpp $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)
	#./$(APP)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

##########################################################
all: exising.out

CC = g++ -g -Wall

exising.out: ising.cpp
	$(CC) -oexising.out ising.cpp -larmadillo -llapack -lblas
