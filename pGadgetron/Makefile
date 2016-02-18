CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC  -std=c++11 -DNDEBUG -I$(IUTILITIES)
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib
#PYTHON = python2.7
#PYINCLUDE = /usr/include/$(PYTHON)
#PYLIBPATH = /lib/$(PYTHON)/config-x86_64-linux-gnu
PYLIB = -l$(PYTHON)

IUTILITIES = ../../iUtilities
LIBIUTIL = ../../build/iUtilities
BUILDPATH = ../../build/xGadgetron

SWIG = swig

LIBS = \
	-lcgt -liutil \
	-lismrmrd \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_system \
	-lpthread

all: pygadgetron

pygadgetron: $(BUILDPATH)/libcgt.a
	$(SWIG) -c++ -python pygadgetron.i
	$(CPP) $(CFLAGS) -c pygadgetron_wrap.cxx -I$(PYINCLUDE)
	$(CPP) $(CFLAGS) -shared pygadgetron_wrap.o -L$(PYLIBPATH) $(PYLIB)\
		 -o $(BUILDPATH)/_pygadgetron.so \
		-L$(BUILDPATH) -L$(LIBIUTIL) $(LIBS) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o

