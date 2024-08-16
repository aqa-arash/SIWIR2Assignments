#CXX := g++ -O3 -Wall -Winline -Wshadow -std=c++17  # Add all the flags here

#OMPFLAGS := -fopenmp
#LDFLAGS := -fopenmp

.PHONY: all clean

all: triangle

triangle: triangle.cpp Colsamm.h
	$(CXX) $(OMPFLAGS) $(LDFLAGS) -o triangle triangle.cpp -lm

clean:
	rm -f triangle
