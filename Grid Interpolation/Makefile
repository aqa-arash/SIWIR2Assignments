CXX := g++ -O3 -Wall -Winline -Wshadow -std=c++17  # Add all the flags here

OMPFLAGS := -fopenmp
LDFLAGS := -fopenmp

.PHONY: all clean

all: intug

intug: intug.cpp VTKparsing.hpp
	$(CXX) $(OMPFLAGS) $(LDFLAGS) -o intug intug.cpp -lm

clean:
	rm -f intug
