CXX := g++ -O3 -Wall -pedantic -std=c++17  # Add all the flags here

OMPFLAGS := -fopenmp
LDFLAGS := -fopenmp

.PHONY: all clean

TARGET := lbm
OBJECTS := $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OMPFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(TARGET).o: $(TARGET).cpp parsing.hpp
	$(CXX) $(OMPFLAGS) -c $<

clean:
	rm -f $(TARGET) $(OBJECTS)
