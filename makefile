CXX = g++
mcg: simulation.h stats.h photon.h vec.h utility.h MCGamma.cpp
	$(CXX) -std=c++17 -static -O3 -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall
debug: simulation.h stats.h photon.h vec.h utility.h MCGamma.cpp
	$(CXX) -std=c++17 -Og -fsanitize=address,undefined -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall
clean:
	rm -f mcg.o
