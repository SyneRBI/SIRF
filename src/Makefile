CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC  -std=c++11 -DNDEBUG -I$(IUTILITIES)
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib
PYTHON = python2.7
PYINCLUDE = /usr/include/$(PYTHON)
PYLIBPATH = /lib/$(PYTHON)/config-x86_64-linux-gnu
PYLIB = -l$(PYTHON)

BUILDPATH = ../../build/xGadgetron
IUTILITIES = ../../iUtilities
LIBIUTIL = ../../build/iUtilities

SWIG = swig

LIBS = \
	-lismrmrd \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_system \
	-lpthread

all: libcgt pygadgetron test_client

libcgt: cgadgetron.o
	$(ARCH) $(ARCHFLAGS) libcgt.a cgadgetron.o
	$(RANLIB) libcgt.a
	mv libcgt.a $(BUILDPATH)

test_client.o: cgadgetron.h xgadgetron.h

test_client: test_client.o libcgt
	$(CPP) $(CFLAGS) -o $(BUILDPATH)/test_client test_client.o \
	-L$(BUILDPATH) -lcgt -L$(LIBIUTIL) -liutil  $(LIBS) $(LDFLAGS)

pygadgetron: cgadgetron.h libcgt
	$(SWIG) -c++ -python pygadgetron.i
	$(CPP) $(CFLAGS) -c pygadgetron_wrap.cxx -I$(PYINCLUDE)
	$(CPP) $(CFLAGS) -shared pygadgetron_wrap.o -L$(PYLIBPATH) $(PYLIB)\
		 -o $(BUILDPATH)/_pygadgetron.so \
		-L$(BUILDPATH) -lcgt -L$(LIBIUTIL) -liutil $(LIBS) $(LDFLAGS)

input:
	ismrmrd_generate_cartesian_shepp_logan -r 2
	mv testdata.h5 ../examples

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f $(BUILDPATH)/test_client

