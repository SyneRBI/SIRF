CFLAGS = -O3 -fPIC -DNDEBUG -DGCC

CGADGETRON = $(SRC_PATH)/xGadgetron/cGadgetron
CGADGETRON_LIBRARY = $(BUILD_PATH)/xGadgetron/libcgt.a

INCLUDE =  -I$(MATLABROOT)/extern/include -I$(CGADGETRON)

LIBS = \
	-lismrmrd \
	-lfftw3 -lfftw3f \
	-lboost_filesystem \
	-lboost_thread \
	-lpthread

TARGET = mgadgetron.mexa64

all: TARGET

OBJECTS = mgadgetron.o

TARGET: $(OBJECTS)
	$(GCC) $(CFLAGS) \
	-shared -Wl,-soname,$(TARGET) \
	-o $(TARGET) $(OBJECTS) \
	$(IUTILITIES_LIBRARY) \
	$(CGADGETRON_LIBRARY) \
	$(LIBS)

%.o: %.c
	$(GCC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.cpp
	$(GCC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	rm -f *.o

