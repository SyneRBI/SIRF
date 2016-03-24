CFLAGS = -O3 -fopenmp -fPIC -DNDEBUG -DGCC

IUTILITIES = $(SRC_PATH)/iUtilities
IUTILITIES_LIBRARY = $(BUILD_PATH)/iUtilities/libiutil.a
CGADGETRON = $(SRC_PATH)/xGadgetron/cGadgetron
CGADGETRON_LIBRARY = $(BUILD_PATH)/xGadgetron/libcgt.a

INCLUDE =  -I$(MATLABROOT)/extern/include \
	 -I$(ISMRMRD_INCLUDE) \
	-I$(IUTILITIES) -I$(CGADGETRON)

LIBS = \
	-lismrmrd \
	-lfftw3 -lfftw3f \
	-lboost_filesystem \
	-lboost_thread \
	-lpthread

TARGET1 = mutilities.mexa64
TARGET2 = mgadgetron.mexa64

all: TARGET1 TARGET2

OBJECTS1 = mutilities.o printer.o
OBJECTS2 = mgadgetron.o

TARGET1: $(OBJECTS1)
	$(GCC) $(CFLAGS) \
	-shared -Wl,-soname,$(TARGET1) \
	-o $(TARGET1) $(OBJECTS1) $(IUTILITIES_LIBRARY)

TARGET2: $(OBJECTS2)
	$(GCC) $(CFLAGS) \
	-shared -Wl,-soname,$(TARGET2) \
	-o $(TARGET2) $(OBJECTS2) \
	$(IUTILITIES_LIBRARY) \
	$(CGADGETRON_LIBRARY) \
	$(LIBS)

%.o: %.c
	$(GCC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.cpp
	$(GCC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	rm -f *.o

