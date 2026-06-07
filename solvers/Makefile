CXX = g++
CXXFLAGS = -O3 -march=native -Wall -std=c++17 -pthread -DDEBUG_RR

SRC = src/main.cpp src/mcsplit.cpp src/mcsplit-dal.cpp src/rrsplit.cpp src/symsplit.cpp src/mcsplit-ll.cpp src/mcsplit-dsb.cpp src/graph.cpp

all: bin/mcsp

bin/mcsp: $(SRC)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o bin/mcsp $(SRC)

clean:
	rm -f bin/mcsp src/*.o