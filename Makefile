# Daniel Graves
# HMACTool

# comment out if not using
#with_mpi = true
#with_boost = true

# set boost include and lib paths
BOOST_ROOT=/usr/local
BOOST_INC_PATH=$(BOOST_ROOT)/include
BOOST_LIB_PATH=$(BOOST_ROOT)/lib

CXX = c++
CXXFLAGS = -O2 --std=c++0x
CPPFLAGS = -I./inc

# mpi
ifdef with_mpi
CXXFLAGS += -DWITH_MPI
MPICXX = mpic++
else
MPICXX = c++
endif

# boost
ifdef with_boost
CXXFLAGS += -DWITH_BOOST
CPPFLAGS += -I$(BOOST_INC_PATH)
LDFLAGS += -L$(BOOST_LIB_PATH) -lboost_filesystem -lboost_system
endif

OBJS = lib/HMAC.opp lib/HMACTool.opp

all: HMACTool

HMACTool: $(OBJS)
	$(MPICXX) $(OBJS) $(LDFLAGS) -o bin/HMACTool

lib/HMAC.opp: inc/HMAC.hpp src/HMAC.cpp
	$(CXX) -c src/HMAC.cpp $(CPPFLAGS) $(CXXFLAGS) -o lib/HMAC.opp

lib/HMACTool.opp: src/HMACTool.cpp
	$(MPICXX) -c src/HMACTool.cpp $(CPPFLAGS) $(CXXFLAGS) -o lib/HMACTool.opp

clean:
	rm -f bin/* lib/*

