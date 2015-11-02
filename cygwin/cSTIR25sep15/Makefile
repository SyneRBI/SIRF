SWIG = /lib/swig-3.0.5/preinst-swig
PYINCLUDE = /usr/include/python2.7
PYLIBPATH = /lib/python2.7/config
PYLIB = -lpython2.7

STIR_IPATH = C:/cygwin64/home/wps46139/Codes/STIR/stirdev-stir-97723408086a/src/include
CSTIR_PATH = C:/Users/wps46139/Codes/PETMR/STIR/cSTIR

CXX = g++
CFLAGS = -O3 -fpermissive -std=c++11 -I$(STIR_IPATH) -I$(CSTIR_PATH)
#CFLAGS = -g -fpermissive -std=c++11 -I$(STIR_IPATH) -I$(CSTIR_PATH)
LDFLAGS = -lrt -lX11

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib
CSTIRLIB = libcstir.a

CSTIR_LIB = -lcstir
#STIR_LIBPATH = /usr/local/lib
STIR_LIBPATH = /home/wps46139/Codes/STIR/Release/lib
#STIR_LIBPATH = /home/wps46139/Codes/STIR/Debug/lib
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

REG_PATH = ./Registries

REG_OBJECTS = \
	$(REG_PATH)/buildblock_registries.o \
	$(REG_PATH)/data_buildblock_registries.o \
	$(REG_PATH)/IO_registries.o \
	$(REG_PATH)/modelling_registries.o \
	$(REG_PATH)/recon_buildblock_registries.o \
	$(REG_PATH)/Shape_buildblock_registries.o \
	$(REG_PATH)/spatial_transformation_registries.o

OBJECTS = ci_stir.o stir_p.o ci_dh.o ci_ex.o ci_tw.o

all: libcstir pystir

libcstir: $(OBJECTS)
	$(ARCH) $(ARCHFLAGS) $(CSTIRLIB) $(OBJECTS)
	$(RANLIB) $(CSTIRLIB)

ci_stir.o: $(CSTIR_PATH)\ci_stir.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

stir_p.o: $(CSTIR_PATH)\stir_p.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_dh.o: $(CSTIR_PATH)\ci_dh.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_ex.o: $(CSTIR_PATH)\ci_ex.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

ci_tw.o: $(CSTIR_PATH)\ci_tw.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

TEST_OBJECTS = main.o cstir_tests.o

main.exe: $(TEST_OBJECTS)
	$(CXX) $(CFLAGS) -o main.exe $(TEST_OBJECTS) $(REG_OBJECTS)\
		-L./ $(CSTIR_LIB) -L$(STIR_LIBPATH)  $(STIR_LIBS) $(LDFLAGS)

pystir: $(CSTIR_PATH)/cstir.h
	$(SWIG) -c++ -python pystir.i
	g++ -c pystir_wrap.cxx -I$(PYINCLUDE) -I$(PYLIBPATH)
	g++ -shared pystir_wrap.o $(REG_OBJECTS) -L$(PYLIBPATH) $(PYLIB) -o _pystir.dll \
		-L./ $(CSTIR_LIB) -L$(STIR_LIBPATH) $(STIR_LIBS) $(LDFLAGS)

%.o: %.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f *.o *.exe
