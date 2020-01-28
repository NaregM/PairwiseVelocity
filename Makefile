# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -I.

# ****************************************************
# Targets needed to bring the executable up to date

main: main.o Cluster.o
	$(CXX) $(CXXFLAGS) -o driver main.o Cluster.o

# The main.o target can be written more simply

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

Cluster.o: Cluster.cpp
	$(CXX) $(CXXFLAGS) -c Cluster.cpp


clean:
	rm -f *~ *.o driver
