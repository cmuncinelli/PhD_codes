all: programs

######################################################################
# Makefile for building PythiaGen - Adapted from the example
######################################################################

# Or just run: g++ PythiaGen_SOFT.cxx -o PythiaGen_SOFT -I/home/daviddc/pythia8/include -O2 -std=c++98 -pedantic -W -Wall -Wshadow -fPIC -L/home/daviddc/pythia8/lib -Wl,-rpath,/home/daviddc/pythia8/lib -lpythia8 -ldl `root-config --cflags` `root-config --libs`
# g++ PythiaGen_SOFT.cxx -o PythiaGen_SOFT -I/home/daviddc/pythia8/include -O2 -std=c++98 -pedantic -W -Wall -Wshadow -fPIC -L/home/daviddc/pythia8/lib -Wl,-rpath,/home/daviddc/pythia8/lib -lpythia8 -ldl `root-config --cflags` `root-config --libs`; g++ PythiaGen_HARD.cxx -o PythiaGen_HARD -I/home/daviddc/pythia8/include -O2 -std=c++98 -pedantic -W -Wall -Wshadow -fPIC -L/home/daviddc/pythia8/lib -Wl,-rpath,/home/daviddc/pythia8/lib -lpythia8 -ldl `root-config --cflags` `root-config --libs`

# Need this to get SHAREDSUFFIX (e.g. dylib or so)
-include ../config.mk

# the compiler to use
CXXCOMP      := g++

PROGS     = PythiaGen
programs: $(PROGS)

# A few variables used in this Makefile:
PYTHIA8      := /home/daviddc/pythia8
STATICLIB    := $(PYTHIA8)/lib/libpythia8.a
SHAREDLIB    := $(PYTHIA8)/lib/libpythia8.$(SHAREDSUFFIX)
DICTCXXFLAGS := -I$(PYTHIA8)/include
ROOTCXXFLAGS := $(DICTCXXFLAGS) $(shell root-config --cflags)
LHAPDF       := /opt/LHAPDF/lib
INTELCXXFLAGS:= -O2 -std=c++98 -pedantic -W -Wall -Wshadow -fPIC -lpythia8 -ldl `root-config --cflags` `root-config --libs`


# Libraries to include if GZIP support is enabled
ifeq (x$(ENABLEGZIP),xyes)
LIBGZIP=-L$(BOOSTLIBLOCATION) -lboost_iostreams -L$(ZLIBLOCATION) -lz
endif

# LDFLAGS1 for static library
LDFLAGS1 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib -lpythia8 -lEG -ldl $(LIBGZIP)


# parnassus
PythiaGen: $(STATICLIB) PythiaGen
	$(CXXCOMP) $(INTELCXXFLAGS) $(ROOTCXXFLAGS) $@.cxx -o $@.exe $(LDFLAGS1)

# Clean up
clean:
	rm -f  *.exe 