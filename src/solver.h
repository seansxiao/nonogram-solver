#include "puzzle.h"
#include <algorithm> 
struct run {
	int s;
	int e;
	int l;
};

// struct cell {
// 	run** runs;
// 	int count;
// 	bool* present;
// };

struct solver {
	int width;
	int height;

	run** row_runs;
	int* row_sizes;

	run** col_runs;
	int* col_sizes;

	// cell** cells;
};

using Solver = solver*;

Solver initialize_solver(Puzzle p);
void initialize_runs(Puzzle p, Solver solv);
void free_solver(Solver solv);

struct state {
	Solver solv;
	Solution solu;
	int row;
	int run_index;
};

using State = state*;

State create_state(Solver s, Solution sol);
void free_state(State st);



bool solve(Puzzle p, Solution s);
bool solve_helper(Puzzle p, State st);
bool solved(Solution s);
