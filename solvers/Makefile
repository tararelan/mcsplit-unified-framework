CXX = g++
CXXFLAGS = -O3 -march=native -Wall -std=c++17 -pthread -DDEBUG_RR

BOOST_VERSION = 1.84.0
BOOST_DIR = boost_1_84_0
BOOST_URL = https://archives.boost.io/release/$(BOOST_VERSION)/source/$(BOOST_DIR).tar.gz
BOOST_INC = $(HOME)/work/$(BOOST_DIR)

SRC = src/main.cpp src/mcsplit.cpp src/mcsplit-dal.cpp src/rrsplit.cpp src/symsplit.cpp src/mcsplit-ll.cpp src/mcsplit-dsb.cpp src/mcsplit-rl.cpp src/graph.cpp src/proof.cpp

all: bin/mcsp

# Downloads and extracts Boost headers if not already present
$(BOOST_INC)/boost/functional/hash.hpp:
	mkdir -p $(HOME)/work
	wget -q -O $(HOME)/work/$(BOOST_DIR).tar.gz $(BOOST_URL)
	tar xzf $(HOME)/work/$(BOOST_DIR).tar.gz -C $(HOME)/work
	rm -f $(HOME)/work/$(BOOST_DIR).tar.gz

bin/mcsp: $(SRC) $(BOOST_INC)/boost/functional/hash.hpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -I $(BOOST_INC) -o bin/mcsp $(SRC)

# Checks the proofs of every algorithm with VeriPB (and with CakePB, if
# CAKE_PB_MCIS is set to the path of cake_pb_mcis)
check-proofs: bin/mcsp
	python3 test/check_proofs.py --solver bin/mcsp $(if $(CAKE_PB_MCIS),--cake $(CAKE_PB_MCIS))

# Checks that mutated solvers and proof loggers are all caught
check-mutants: $(BOOST_INC)/boost/functional/hash.hpp
	CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) -I $(BOOST_INC)" python3 test/mutants.py

clean:
	rm -f bin/mcsp src/*.o

.PHONY: all check-proofs check-mutants clean