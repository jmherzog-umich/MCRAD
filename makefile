CXX = g++
mcg:
	$(CXX) -std=c++14 -static -o mcg MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h -Wall
clean:
	rm -f mcg
