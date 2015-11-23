CC=mpiCC
CXXFLAGS= -O3
LDFLAGS=
SOURCES=apspmpi.cpp apspmpi_driver.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mpiout

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CXXFLAGS) $< -o $@

clean:
	rm mpiout 
