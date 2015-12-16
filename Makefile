CC=mpiCC
CXXFLAGS= -O3
LDFLAGS=
SOURCES=apspmpi.cpp apspmpi_solver.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=apspmpi

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.o:
	$(CC) $(CXXFLAGS) $< -o $@

clean:
	rm apspmpi
