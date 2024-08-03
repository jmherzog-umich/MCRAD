CXX = g++
CXXW32 = i686-w64-mingw32-c++
CXXW64 = x86_64-w64-mingw32-c++
CPPFLAGS = 
MPICPPFLAGS = -D_MPI -L /usr/lib64/openmpi/lib
MPICXX = mpicxx
SOURCES = $(wildcard *.h) $(wildcard **/*.h)

mcrad: mcrad.cpp $(SOURCES)
	$(CXX) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp $(SOURCES) -Wall
	
win32: mcrad.cpp $(SOURCES)
	$(CXXW32) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.W32.exe mcrad.cpp $(SOURCES) -Wall
	
win64: mcrad.cpp $(SOURCES)
	$(CXXW64) $(CPPFLAGS) -std=c++17 -static -Ofast -o mcrad.W64.exe mcrad.cpp $(SOURCES) -Wall
	
mpi: mcrad.cpp $(SOURCES)
	$(MPICXX) $(MPICPPFLAGS) -std=c++17 -static -Ofast -o mcrad.o mcrad.cpp $(SOURCES) -Wall
	
debug: mcrad.cpp $(SOURCES)
	$(CXX) $(CPPFLAGS) -std=c++17 -Og -fsanitize=address,undefined -o mcrad.o mcrad.cpp $(SOURCES) -Wall -Wuninitialized -Wextra -g -rdynamic -D_DEBUG
	
clean:
	rm -f mcrad.o mcrad.W32.exe mcrad.W64.exe
