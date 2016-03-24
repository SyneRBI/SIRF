CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

LIBIUTIL = $(BUILD_PATH)/iUtilities

all: libiutil

libiutil: ci_tw.o data_handle.o text_writer.o
	$(ARCH) $(ARCHFLAGS) libiutil.a ci_tw.o data_handle.o text_writer.o
	$(RANLIB) libiutil.a
	mv libiutil.a $(LIBIUTIL)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
