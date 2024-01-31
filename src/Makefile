

CXX = g++
CXXFLAGS = -c -std=c++0x -O3 -w

TEST = Test.o
OBJ = Vertex.o Graph.o Matrix.o MinHeap.o SCGTree.o AuxiliaryData.o


test: $(OBJ) $(TEST)
	$(CXX) $(OBJ) $(TEST) $(INCLUDE) $(LIBMETIS) -o test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBMETIS) $< -o $@




.PHONY: clean
clean:
	rm -f *.o test BAY.scg