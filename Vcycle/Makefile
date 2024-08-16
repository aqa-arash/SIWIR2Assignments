CXX := g++ -O3 -Wall -Winline -Wshadow -std=c++17  # Add all the flags here

OMPFLAGS := -fopenmp
LDFLAGS := -fopenmp

.PHONY: all clean

all: mgsolve

mgsolve: mgsolve.cpp gssolve.hpp
	$(CXX) $(OMPFLAGS) $(LDFLAGS) -o mgsolve mgsolve.cpp -lm

clean:
	rm -f mgsolve
