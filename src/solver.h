#include "puzzle.h"

struct run {
	int s;
	int e;
	int l;

	int start(int sNew) {
		if (e - sNew + 1 < l) {
			return CONFLICT;
		}
		else if (sNew > s) {
			s = sNew;
			return PROGRESS;
		}
		else {
			return SAME;
		}
	}

	int end(int eNew) {
		if (eNew - s + 1 < l) {
			return CONFLICT;
		}
		else if (eNew < e) {
			e = eNew;
			return PROGRESS;
		}
		else {
			return SAME;
		}
	}
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
bool filled(Solution s);
