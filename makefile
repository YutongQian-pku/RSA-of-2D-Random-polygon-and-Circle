main: main.cpp particle.o voxel.o cell.o 
	g++ -std=c++17 -O3 main.cpp particle.o voxel.o cell.o  -o main
particle.o: particle.cpp particle_type.h
	g++ -std=c++17 -O3 -c particle.cpp
voxel.o: voxel.cpp
	g++ -std=c++17 -O3 -c voxel.cpp
cell.o: cell.cpp
	g++ -std=c++17 -O3 -c cell.cpp
clean:
	rm *.o main