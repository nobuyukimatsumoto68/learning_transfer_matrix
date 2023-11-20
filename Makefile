CXX      = g++
CXXFLAGS = -O3 -std=c++14 -lfftw3 -fopenmp
INCLUDES = -I/usr/local/include
INCLUDES += -L/usr/local/lib/


a.out:	main.cpp header.h
	$(CXX) main.cpp $(CXXFLAGS) $(INCLUDES)
