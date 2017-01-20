CPP = /usr/bin/g++

all : c2m.exe

OBJECTS = main.o c2m.o

c2m.exe: $(OBJECTS)
	$(CPP) -o c2m.exe $(OBJECTS)

%.o: %.cpp
	$(CPP) -c -o $@ $<

clean:
	rm -f *.o

