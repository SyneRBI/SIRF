CPP = /usr/bin/g++
#CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(ISMRMRD_HOME)/include
CFLAGS = -O3 -fPIC -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(ISMRMRD_HOME)/include
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = $(SRC_PATH)/iUtilities

all: libcgt

libcgt: cgadgetron.o ismrmrd_fftw.o noise_filter.o
	$(ARCH) $(ARCHFLAGS) libcgt.a cgadgetron.o ismrmrd_fftw.o noise_filter.o
	$(RANLIB) libcgt.a
	mv libcgt.a $(BUILD_PATH)/xGadgetron

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
