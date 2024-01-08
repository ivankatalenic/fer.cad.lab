CXX=g++
CXXFLAGS=-std=c++17 -Wall -O2
CPPFLAGS=
LDFLAGS=
LDLIBS=

OBJECTS := $(patsubst %.cc,%.o,$(wildcard *.cc))

.PHONY: build
build: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o program $(OBJECTS)

.PHONY: run
run: build
	./program.exe
