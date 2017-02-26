.PHONY: cleanall, mrproper

.SUFFIXES:

# Compile parameters
COMP=g++
FLAGS= -std=c++0x -W -Wall -g --pedantic -DNDEBUG

all: sph

# Objects used
sph: SPH.o Geometry.o ParaView.o Neighborhood.o #Neighborhood_performance.o
	$(COMP) $^ $(FLAGS) -o $@

SPH.o: SPH.hpp #SPH.cpp Geometry.cpp ParaView.cpp

#Neighborhood_performance.o: SPH.hpp

Neighborhood.o: SPH.hpp

Geometry.o: SPH.hpp

ParaView.o: SPH.hpp

# Compilation :
%.o: %.cpp
	$(COMP) -c $< $(FLAGS) -o $@

cleanall:
	-rm -rf *.o

mrproper: cleanall
	-rm -rf sph
# End of Makefile
