CPP = g++
CFLAGS = -O3 -fPIC -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(ISMRMRD_HOME)/include
#LDFLAGS = -lrt
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib
PYLIB = -l$(PYTHON)

IUTILITIES = $(SRC_PATH)/iUtilities
LIBIUTIL = $(BUILD_PATH)/iUtilities
BUILDPATH = $(BUILD_PATH)/xGadgetron
CGADGETRON = ../cGadgetron

SWIG = swig

LIBS = \
	-lcgt -liutil \
	-lismrmrd \
	-lfftw3 -lfftw3f \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_system \
	-lpthread

all: pygadgetron

pygadgetron: $(BUILDPATH)/libcgt.a
	$(SWIG) -c++ -python -I$(CGADGETRON) pygadgetron.i
	$(CPP) $(CFLAGS) -c pygadgetron_wrap.cxx -I$(PYINCLUDE) -I$(CGADGETRON)
	$(CPP) $(CFLAGS) -shared pygadgetron_wrap.o -L$(PYLIBPATH) $(PYLIB)\
		 -o $(BUILDPATH)/_pygadgetron.so \
		-L$(BUILDPATH) -L$(LIBIUTIL) $(LIBS) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o

