CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(ISMRMRD_HOME)/include
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = ../../iUtilities
BUILDPATH = ../../build/xGadgetron

all: libcgt

libcgt: cgadgetron.o ismrmrd_fftw.o
	$(ARCH) $(ARCHFLAGS) libcgt.a cgadgetron.o ismrmrd_fftw.o
	$(RANLIB) libcgt.a
	mv libcgt.a $(BUILDPATH)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
