CXX = g++
CFLAGS = -O3 -fpermissive -fPIC -std=c++11 -I$(STIR_SRC)/include -I$(CSTIR_SRC)
LDFLAGS = $(X11LIB) -lrt

SWIG = swig

ifneq ($(wildcard /*/lib/X11/.*),)
	X11LIB = -lX11
endif

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

CSTIRLIB = libcstir.a
CSTIR_LIB = -lcstir

STIR_LIBS = \
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

OBJECTS = ci_stir.o stir_p.o ci_dh.o ci_ex.o ci_tw.o

all: libcstir pystir

libcstir: $(OBJECTS)
	$(ARCH) $(ARCHFLAGS) $(CSTIRLIB) $(OBJECTS)
	$(RANLIB) $(CSTIRLIB)

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

ci_stir.o: $(CSTIR_SRC)/ci_stir.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

stir_p.o: $(CSTIR_SRC)/stir_p.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_dh.o: $(CSTIR_SRC)/ci_dh.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_ex.o: $(CSTIR_SRC)/ci_ex.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_tw.o: $(CSTIR_SRC)/ci_tw.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

TEST_OBJECTS = main.o cstir_tests.o

main.exe: $(TEST_OBJECTS)
	$(CXX) $(CFLAGS) -o main.exe $(TEST_OBJECTS) $(REG_OBJECTS)\
		-L./ $(CSTIR_LIB) $(STIR_LIBS) $(LDFLAGS)

pystir: $(CSTIR_SRC)/cstir.h $(REG_OBJECTS) $(CSTIRLIB)
	$(SWIG) -c++ -python pystir.i
	$(CXX) $(CFLAGS) -c pystir_wrap.cxx -I$(PYINCLUDE) -I$(PYLIBPATH)
	$(CXX) $(CFLAGS) -shared pystir_wrap.o $(REG_OBJECTS) -L$(PYLIBPATH) -l$(PYTHON) -o _pystir.so \
		-L./ $(CSTIR_LIB) -L$(STIR_LIBPATH) $(STIR_LIBS) $(LDFLAGS)

%.o: %.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.exe
