all: default ref

default: main.cpp
	g++ -std=c++11 -fopenmp -O3 -g -o nonogram_solver main.cpp puzzle.cpp para.cpp

ref: main.cpp
	g++ -std=c++11 -fopenmp -O3 -g -o nonogram_solver_ref main.cpp puzzle.cpp solver.cpp
