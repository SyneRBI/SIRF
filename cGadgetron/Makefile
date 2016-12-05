CPP = /usr/bin/g++
CFLAGS = -O3 -fPIC -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(ISMRMRD_HOME)/include -I$(HDF5_INCLUDE)
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = $(SRC_PATH)/iUtilities

all: libcgt

libcgt: cgadgetron.o xgadgetron.o ismrmrd_fftw.o
	$(ARCH) $(ARCHFLAGS) libcgt.a cgadgetron.o xgadgetron.o ismrmrd_fftw.o
	$(RANLIB) libcgt.a
	mv libcgt.a $(BUILD_PATH)/xGadgetron
	
HEADERS = xgadgetron.h gadgetron_data_containers.h

%.o: %.cpp $(HEADERS) $(BUILD_PATH)/iUtilities/libiutil.a
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
