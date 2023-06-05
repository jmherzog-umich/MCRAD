SET PATH=C:\Users\jmherzog\Software\mingw64\bin\;%PATH%
g++ -std=c++14 -static -o mcg.exe MCGamma.cpp photon.h stats.h vec.h simulation.h utility.h -Wall
pause
