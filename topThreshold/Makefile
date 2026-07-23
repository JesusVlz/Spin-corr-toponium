# Makefile for main_topThreshold_lhe.cc
# Standalone build, independent of the Pythia8 examples/ directory.
#
# Two ways to point this at your Pythia8 installation:
#
# 1) If pythia8-config is available (installed via "make install" in the
#    Pythia8 build, ends up in PREFIX_BIN/bin), it is used automatically:
#      make
#
# 2) Otherwise, set PYTHIA8 to the Pythia8 install/build directory
#    (the one containing include/ and lib/), e.g. the top of the Pythia8
#    source tree if you never ran "make install":
#      make PYTHIA8=/path/to/pythia8
#
# Clean:
#      make clean

# EDIT THIS: path to your Pythia8 build/install directory (the one that
# contains include/Pythia8/Pythia.h and lib/libpythia8.*). Fill in manually.
PYTHIA8 := /home/jesusavc/hep-software/pythia8317/

CXX        := g++
CXX_COMMON := -O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC

EXE := main_topThreshold_lhe
SRC := main_topThreshold_lhe.cc

PYTHIA8_CONFIG := $(shell command -v pythia8-config 2> /dev/null)

ifneq ($(PYTHIA8_CONFIG),)
  PYTHIA8_CXXFLAGS := $(shell pythia8-config --cppflags)
  PYTHIA8_LIBS     := $(shell pythia8-config --libs)
else
  PYTHIA8_CXXFLAGS := -I$(PYTHIA8)/include
  PYTHIA8_LIBS     := -L$(PYTHIA8)/lib -Wl,-rpath,$(PYTHIA8)/lib -lpythia8 -ldl
endif

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC)
	$(CXX) $(CXX_COMMON) $(PYTHIA8_CXXFLAGS) $(SRC) -o $(EXE) $(PYTHIA8_LIBS)

clean:
	rm -f $(EXE) *.lhe