CC=`root-config --cxx`
CFLAGS=-c -g -Wall `root-config --cflags`
LDFLAGS=-lTreePlayer `root-config --ldflags`
LDLIBS=-lTreePlayer `root-config --glibs`
SOURCES=HitValidation.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=HitValidation

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ./*~ ./*.o ./HitValidation
