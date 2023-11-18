CXX      = g++
CXXFLAGS = -O2 -std=c++11
INCLUDES = -I/usr/local/include


a.out:	main.cpp header.h
	$(CXX) main.cpp $(CXXFLAGS) $(INCLUDES)
