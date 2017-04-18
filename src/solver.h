#include "puzzle.h"

struct run {
	int s;
	int e;
	int l;
};

struct cell {
	run** runs;
	int count;
	bool* present;
};

struct solver {
	int width;
	int height;

	run** row_runs;
	int* row_sizes;

	run** col_runs;
	int* col_sizes;

	cell** cells;
};

using Solver = solver*;

struct state {
	Solver sol;
	Solution s;
};

using State = state*;



Solver initialize_solver(Puzzle p);
void initialize_runs(Puzzle p, Solver sol);
State create_state(Solver s, Solution sol);
void solve(Puzzle p, Solution s);
