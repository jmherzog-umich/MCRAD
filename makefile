CXX = g++
mcg: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall
debug: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) -std=c++17 -Og -fsanitize=address,undefined -o mcrad.o mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h -Wall -g -rdynamic -D_DEBUG
clean:
	rm -f mcrad.o
