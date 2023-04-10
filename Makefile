# Simple Makefile

CXX = g++
CXXFLAGS = -Wall -ggdb `root-config --cflags` `fastjet-config --cxxflags` 
LDFLAGS = `root-config --libs` `fastjet-config --libs --plugins` 

all: runComp.exe 

runComp.exe: jetComparison.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

jetComparison.o: jetComparison.cxx 

clean:
	$(RM) *.o *~ *.exe

%.o : %.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $<


