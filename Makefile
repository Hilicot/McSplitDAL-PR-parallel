CXX := g++
CXXFLAGS := -O2 -march=native -fopenmp
all: dal

prelim:
	mkdir -p ./build
	
dal: prelim mcsplitDAL+PR-iter.cpp graph.cpp graph.h mcs.h mcs.cpp stats.h args.h test_utility.cpp reward.cpp reward.h $(shell find heuristics -type f)
	$(CXX) $(CXXFLAGS) -Wall -std=c++2a -o build/mcsplit-dal mcsplitDAL+PR-iter.cpp graph.cpp mcs.h mcs.cpp test_utility.cpp reward.cpp $(shell find heuristics -type f -name '*.cpp') -pthread

debug: prelim mcsplitDAL+PR-iter.cpp graph.cpp graph.h mcs.h mcs.cpp stats.h args.h test_utility.cpp reward.cpp reward.h $(shell find heuristics -type f)
	$(CXX) $(CXXFLAGS) -pg -g -Wall -std=c++2a -o build/mcsplit-dal mcsplitDAL+PR-iter.cpp graph.cpp mcs.h mcs.cpp test_utility.cpp reward.cpp $(shell find heuristics -type f -name '*.cpp') -pthread

leak: prelim mcsplitDAL+PR-iter.cpp graph.cpp graph.h mcs.h mcs.cpp stats.h args.h test_utility.cpp reward.cpp reward.h $(shell find heuristics -type f)
	$(CXX) $(CXXFLAGS) -g -fsanitize=address -Wall -std=c++2a -o build/mcsplit-dal mcsplitDAL+PR-iter.cpp graph.cpp mcs.h mcs.cpp test_utility.cpp reward.cpp $(shell find heuristics -type f -name '*.cpp') -pthread

clean:
	rm -rf build
