CXX = g++
CXXFLAGS = -O3 -march=native -Wall -std=c++17 -pthread

SRC = src/main.cpp src/mcsplit.cpp src/mcsplit-dal.cpp src/graph.cpp

all: bin/mcsp

bin/mcsp: $(SRC)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o bin/mcsp $(SRC)

clean:
	rm -f bin/mcsp src/*.o