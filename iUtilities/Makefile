CPP = g++
CFLAGS = -O3 -fopenmp -fPIC -std=c++11 -DNDEBUG
LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

SWIG = swig

IUTILITIES = $(SRC_PATH)/iUtilities
LIBIUTIL = $(BUILD_PATH)/iUtilities
INCLUDE =  -I$(MATLABROOT)/extern/include -I$(IUTILITIES)

ifneq ($(MATLABROOT),)
all: libiutil.a _pyiutil.so mutilities.mexa64
else
all: libiutil.a _pyiutil.so
endif

libiutil.a: iutilities.o
	$(ARCH) $(ARCHFLAGS) libiutil.a iutilities.o
	$(RANLIB) libiutil.a
	mv libiutil.a $(LIBIUTIL)

_pyiutil.so: $(LIBIUTIL)/libiutil.a
	$(SWIG) -c++ -python pyiutil.i
	$(CPP) $(CFLAGS) -c pyiutil_wrap.cxx -I$(PYINCLUDE)
	$(CPP) $(CFLAGS) -shared  -o _pyiutil.so pyiutil_wrap.o -L$(PYLIBPATH) -l$(PYTHON) \
		$(LIBIUTIL)/libiutil.a \
		$(LDFLAGS)

ifneq ($(MATLABROOT),)
mutilities.mexa64: mutilities.o $(LIBIUTIL)/libiutil.a
	$(GCC) $(CFLAGS) \
	-shared -Wl,-soname,mutilities.mexa64 \
	-o mutilities.mexa64 mutilities.o $(LIBIUTIL)/libiutil.a
endif

%.o: %.c
	$(GCC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.cpp data_handle.h
	$(CPP) $(CFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.cxx
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
