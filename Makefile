CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

LIBIUTIL = $(BUILD_PATH)/iUtilities

all: libiutil

libiutil: data_handle.o
	$(ARCH) $(ARCHFLAGS) libiutil.a data_handle.o
	$(RANLIB) libiutil.a
	mv libiutil.a $(LIBIUTIL)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
