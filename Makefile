CXX := g++
CXXFLAGS := -O2 -march=native -fopenmp
all: dal

prelim:
	mkdir -p ./build

rec: prelim mcsp_rec.cpp graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++2a -o build/recur graph.cpp mcsp_rec.cpp test_utility.cpp -pthread

iter: prelim mcsp_iter.cpp graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++2a -o build/iter graph.cpp mcsp_iter.cpp test_utility.cpp -pthread
	
dal: prelim mcsplit+DAL.cpp graph.cpp graph.h mcs.h mcs.cpp stats.h args.h test_utility.cpp reward.cpp reward.h $(shell find heuristics -type f)
	$(CXX) $(CXXFLAGS) -Wall -std=c++2a -o build/mcsplit-dal mcsplit+DAL.cpp graph.cpp mcs.h mcs.cpp test_utility.cpp reward.cpp $(shell find heuristics -type f -name '*.cpp') -pthread

clean:
	rm -rf build
