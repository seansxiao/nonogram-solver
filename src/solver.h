#include "puzzle.h"

struct run {
	int s;
	int e;
	int l;

	inline int start(int sNew) {
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

	inline int end(int eNew) {
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

struct solver {
	int width;
	int height;

	run** row_runs;
	int* row_sizes;
	bool* solved_rows;

	run** col_runs;
	int* col_sizes;
	bool* solved_cols;
};

using Solver = solver*;

inline Solver initialize_solver(Puzzle p);
inline void initialize_runs(Puzzle p, Solver solv);
inline void free_solver(Solver solv);

struct state {
	Solver solv;
	Solution solu;
	int row;
	int run_index;
};

using State = state*;

inline State create_state(Solver s, Solution sol);
inline void free_state(State st);



bool solve(Puzzle p, Solution s);
bool solve_helper(Puzzle p, State st);
inline bool filled(Solution s);
