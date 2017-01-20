CXX = g++
EXTRA_CFLAGS=-O3 -DNDEBUG
#EXTRA_CFLAGS=-g
CFLAGS = $(EXTRA_CFLAGS) -fpermissive -fPIC -std=c++11 -I$(IUTILITIES) -I$(STIR_SRC)/include -I$(CSTIR_SRC)

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

IUTILITIES = $(SRC_PATH)/iUtilities

CSTIRLIB = libcstir.a
CSTIR_LIB = -lcstir

OBJECTS = cstir.o cstir_p.o cstir_tw.o

all: libcstir.a

libcstir.a: $(OBJECTS)
	$(ARCH) $(ARCHFLAGS) $(CSTIRLIB) $(OBJECTS)
	$(RANLIB) $(CSTIRLIB)
	mv libcstir.a $(BUILD_PATH)/xSTIR

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a
