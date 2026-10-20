.PHONY: all clean

CXX:=g++
CXXFLAGS:=-O3 -std=c++11 -fopenmp -g

vg:=../bin/nvg

all: test clean

test: build_graph $(vg)
	prove -v t

$(vg):
	cd .. && $(MAKE) bin/nvg

build_graph: build_graph.cpp
	cd .. && $(MAKE) test/build_graph

clean:
	rm -f build_graph
