CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
TARGET = tsp_solver

SOURCES = main.cc mst.cc karp.cc mine.cc
OBJECTS = $(SOURCES:.cc=.o)
DATASETS = a280.tsp xql662.tsp kz9976.tsp mona-lisa100K.tsp

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET) tsp_results.csv

test-all: $(TARGET)
	@if [ ! -d "data" ]; then echo "Error: 'data' not found!"; exit 1; fi
	./$(TARGET) $(patsubst %,data/%,$(DATASETS))

benchmark: test-all

.PHONY: all clean test-all benchmark
