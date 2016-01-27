CPP = /usr/bin/g++
CFLAGS = -O3 -fopenmp -fPIC  -std=c++11 -DNDEBUG -I$(IUTILITIES) -I$(CGADGETRON)
LDFLAGS = -lrt

BUILDPATH = $(XGADGETRON_BUILD)

LIBS = \
	-lcgt -liutil \
	-lismrmrd \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_system \
	-lpthread

all: test_client

test_client: test_client.o
	$(CPP) $(CFLAGS) -o $(BUILDPATH)/test_client test_client.o \
	-L$(BUILDPATH) -L$(LIBIUTIL) $(LIBS) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f $(BUILDPATH)/test_client
