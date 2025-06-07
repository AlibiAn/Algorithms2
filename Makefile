CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
DEBUG_FLAGS = -std=c++17 -Wall -Wextra -g

TARGET = tsp_solver
SOURCES = main.cpp mst.cc held_karp.cc results_tracker.cc
HEADERS = mst.h held_karp.h results_tracker.h
OBJECTS = $(SOURCES:.cpp=.o) $(SOURCES:.cc=.o)
DATASETS = a280.tsp xql662.tsp kz9976.tsp mona-lisa100K.tsp

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

debug: CXXFLAGS = $(DEBUG_FLAGS)
debug: clean $(TARGET)

clean:
	rm -f *.o $(TARGET) results.csv results_*.csv

test-all: $(TARGET)
	@echo "Dataset,Algorithm,Vertices,Tour_Cost,Runtime_ms,Approximation_Ratio" > results.csv
	@for dataset in $(DATASETS); do \
		if [ -f "data/$$dataset" ]; then \
			./$(TARGET) data/$$dataset --csv >> results.csv; \
		fi; \
	done

benchmark: test-all

setup:
	@mkdir -p data

submit: clean
	@tar -czf cse331_tsp_assignment.tar.gz *.cpp *.cc *.h Makefile

verify: $(TARGET)
	@if [ -f "data/a280.tsp" ]; then \
		echo "Dataset,Algorithm,Vertices,Tour_Cost,Runtime_ms,Approximation_Ratio" > results.csv; \
		./$(TARGET) data/a280.tsp --csv >> results.csv; \
	fi

.PRECIOUS: %.o
.PHONY: all debug clean test-all benchmark setup submit verify
