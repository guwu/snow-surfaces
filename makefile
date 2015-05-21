
all: Main.cpp MarchingCubes.cpp DataStructs.cpp Reader.cpp
	g++ -o snow_surface -O3 -lm -fopenmp Main.cpp MarchingCubes.cpp DataStructs.cpp Reader.cpp


clean:
	rm -f snow_surface