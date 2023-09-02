CXX = g++
mcg: MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) -std=c++17 -static -Ofast -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall
debug: MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) -std=c++17 -Og -fsanitize=address,undefined -o mcg.o MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall -g -rdynamic -D_DEBUG
clean:
	rm -f mcg.o
