CC=mpiCC
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=apsp-mpi.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=mpiout

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm mpiout 
