CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

INCLUDE =  -I$(MATLABROOT)/extern/include -I$(IUTILITIES)
LIBIUTIL = $(BUILD_PATH)/iUtilities

all: libiutil mutilities.mexa64

libiutil: iutilities.o
	$(ARCH) $(ARCHFLAGS) libiutil.a data_handle.o
	$(RANLIB) libiutil.a
	mv libiutil.a $(LIBIUTIL)

mutilities.mexa64: mutilities.o
	$(GCC) $(CFLAGS) \
	-shared -Wl,-soname,mutilities.mexa64 \
	-o mutilities.mexa64 mutilities.o $(LIBIUTIL)/libiutil.a

%.o: %.cpp data_handle.h
	$(CPP) $(CFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	rm -f *.o
