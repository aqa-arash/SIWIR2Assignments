CXX := g++ -O3 -Wall -std=c++17 -pedantic  -fopenmp  # Add all the flags here

.PHONY: all clean

all: waveguide

waveguide: waveguide.cpp ./Source/Colsamm.h inversepower.hpp io.hpp
	$(CXX)-o waveguide waveguide.cpp -lm

clean:
	rm -f waveguide
