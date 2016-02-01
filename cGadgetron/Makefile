CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC  -std=c++11 -DNDEBUG -I$(IUTILITIES)
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = ../../iUtilities
BUILDPATH = ../../build/xGadgetron

LIBS = \
	-lismrmrd \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_system \
	-lpthread

all: libcgt

libcgt: cgadgetron.o
	$(ARCH) $(ARCHFLAGS) libcgt.a cgadgetron.o
	$(RANLIB) libcgt.a
	mv libcgt.a $(BUILDPATH)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
