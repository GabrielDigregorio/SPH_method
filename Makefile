.PHONY: cleanall, mrproper

.SUFFIXES:

# Compile parameters
COMP=g++
FLAGS= -std=c++0x -W -Wall -g --pedantic -DNDEBUG

all: sph

# Objects used
sph: SPH.o Neighborhood_performance.o Geometry.o ParaView.o
	$(COMP) $^ $(FLAGS) -o $@

main.o: SPH.hpp

Neighborhood_performance.o: SPH.hpp

Geometry.o: SPH.hpp

Paraview.o: SPH.hpp

# Compilation :
%.o: %.cpp
	$(COMP) -c $< $(FLAGS) -o $@

cleanall:
	-rm -rf *.o

mrproper: cleanall
	-rm -rf sph
# End of Makefile
