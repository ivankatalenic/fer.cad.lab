CXX=g++
CXXFLAGS=-std=c++17 -Wall
LDFLAGS=
LDLIBS=
LDPATHS=
LINKSFLAGS=

OBJECTS := $(patsubst %.cc,%.o,$(wildcard *.cc))

.PHONY: build
build: $(OBJECTS)
	$(CXX) -o program $(OBJECTS) $(LINKFLAGS) $(LDPATHS) $(LDLIBS)

.PHONY: run
run: build
	./program.exe
