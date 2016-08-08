CXX = g++
EXTRA_CFLAGS=-O3 -DNDEBUG
#EXTRA_CFLAGS=-g
CFLAGS = $(EXTRA_CFLAGS) -fpermissive -fPIC -std=c++11 -I$(IUTILITIES) -I$(STIR_SRC)/include -I$(CSTIR_SRC)
# enable this (or use it on the command line) if you compiled STIR with OPENMP=ON
OPENMPFLAGS=-fopenmp
#EXTRA_LDFLAGS=-g
LDFLAGS = $(X11LIB) -lrt $(OPENMPFLAGS) $(EXTRA_LDFLAGS)

BUILDPATH = $(BUILD_PATH)/xSTIR

SWIG = swig

ifneq ($(wildcard /*/lib/X11/.*),)
	X11LIB = -lX11
endif

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = $(SRC_PATH)/iUtilities
LIBIUTIL = $(BUILD_PATH)/iUtilities
IUTIL_LIB = -liutil

CSTIRLIB = libcstir.a
CSTIR_LIB = -lcstir

STIR_LIBS = \
			-lanalytic_FBP2D \
			-lanalytic_FBP3DRP \
			-literative_OSMAPOSL \
			-literative_OSSPS \
			-lrecon_buildblock \
			-lIO \
			-ldata_buildblock \
			-ldisplay \
			-leval_buildblock \
			-llistmode_buildblock \
			-lmodelling_buildblock \
			-lscatter_buildblock \
			-lShape_buildblock \
			-lnumerics_buildblock \
			-lspatial_transformation_buildblock \
			-lbuildblock

REG_OBJECTS = \
	buildblock_registries.o \
	data_buildblock_registries.o \
	IO_registries.o \
	modelling_registries.o \
	recon_buildblock_registries.o \
	Shape_buildblock_registries.o \
	spatial_transformation_registries.o

all: _pystir.so

buildblock_registries.o: $(STIR_SRC)/buildblock/buildblock_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

data_buildblock_registries.o: $(STIR_SRC)/data_buildblock/data_buildblock_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

IO_registries.o: $(STIR_SRC)/IO/IO_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

modelling_registries.o: $(STIR_SRC)/modelling_buildblock/modelling_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

recon_buildblock_registries.o: $(STIR_SRC)/recon_buildblock/recon_buildblock_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

Shape_buildblock_registries.o: $(STIR_SRC)/Shape_buildblock/Shape_buildblock_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

spatial_transformation_registries.o: \
		$(STIR_SRC)/spatial_transformation_buildblock/spatial_transformation_registries.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

_pystir.so: $(CSTIR_SRC)/cstir.h $(REG_OBJECTS) $(BUILDPATH)/$(CSTIRLIB)
	$(SWIG) -c++ -python pystir.i
	$(CXX) $(CFLAGS) -c pystir_wrap.cxx -I$(PYINCLUDE) -I$(PYLIBPATH)
	$(CXX) $(CFLAGS) -shared  -o _pystir.so pystir_wrap.o $(REG_OBJECTS) -L$(PYLIBPATH) -l$(PYTHON) \
		-L$(LIBIUTIL) $(IUTIL_LIB) \
		-L$(BUILDPATH) $(CSTIR_LIB) -L$(STIR_LIBPATH) $(STIR_LIBS) $(LDFLAGS)

%.o: %.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.exe *.a *.pyc
