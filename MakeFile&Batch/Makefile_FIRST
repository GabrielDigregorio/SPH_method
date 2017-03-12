.PHONY: cleanall, mrproper

.SUFFIXES:

# Compile parameters
COMP=g++
FLAGS= -std=c++11  -W -Wall -g --pedantic

all: sph

# Objects used
sph: Geometry.o ParaView.o Neighborhood.o MemUsage.o Neighborhood_performance.o #SPH.o
	$(COMP) $^ $(FLAGS) -o $@

#SPH.o: SPH.hpp
Neighborhood_performance.o: SPH.hpp

Neighborhood.o: SPH.hpp

Geometry.o: SPH.hpp

ParaView.o: SPH.hpp

MemUsage.o: SPH.hpp

# Compilation :
%.o: %.cpp
	$(COMP) -c $< $(FLAGS) -o $@

cleanall:
	-rm -rf *.o

mrproper: cleanall
	-rm -rf sph
# End of Makefile
