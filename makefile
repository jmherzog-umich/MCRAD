CXX = g++
mcg: simulation.h stats.h photon.h vec.h utility.h MCGamma.cpp
	$(CXX) -std=c++14 -static -O3 -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h -Wall
debug: simulation.h stats.h photon.h vec.h utility.h MCGamma.cpp
	$(CXX) -std=c++14 -static -Og -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h -Wall
clean:
	rm -f mcg.o
