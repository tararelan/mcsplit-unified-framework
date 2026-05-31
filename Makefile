CXX = g++
CXXFLAGS = -O3 -march=native -Wall -std=c++17 -pthread

all: mcsp

mcsp: main.cpp mcsplit.cpp graph.cpp
	$(CXX) $(CXXFLAGS) -o mcsp main.cpp mcsplit.cpp graph.cpp

clean:
	rm -f mcsp *.o