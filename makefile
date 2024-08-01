CXX = g++
CXXW32 = i686-w64-mingw32-c++
CXXW64 = x86_64-w64-mingw32-c++
CPPFLAGS = 
MPICPPFLAGS = -D_MPI -L /usr/lib64/openmpi/lib
MPICXX = mpicxx

mcrad: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h spectrum.h -Wall
	
win32: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXXW32) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.W32.exe mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h spectrum.h -Wall
	
win64: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXXW64) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.W64.exe mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h spectrum.h -Wall
	
mpi: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(MPICXX) $(MPICPPFLAGS) -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h spectrum.h -Wall
	
debug: mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h
	$(CXX) $(CPPFLAGS) -std=c++17 -Og -fsanitize=address,undefined -o mcrad.o mcrad.cpp photon.h stats.h vec.h simulation.h utility.h medium.h beam.h camera.h image.h raypath.h spectrum.h -Wall -Wuninitialized -Wextra -g -rdynamic -D_DEBUG
	
clean:
	rm -f mcrad.o mcrad.W32.exe mcrad.W64.exe
