#include "solver.h"

#include <stdio.h>
#include <stdlib.h>

Solver initialize_solver(Puzzle p) {
	solver* sol = (struct solver*)(malloc(sizeof(struct solver)));

	sol->width = p->width;
	sol->height = p->height;
	
	sol->row_sizes = p->row_sizes;
	sol->col_sizes = p->col_sizes;

	sol->row_runs = (struct run**)(malloc(sizeof(struct run*) * sol->height));
	for (int i = 0; i < p->height; i++) {
		sol->row_runs[i] = (struct run*)(malloc(sizeof(struct run) * sol->row_sizes[i]));
	}

	sol->col_runs = (struct run**)(malloc(sizeof(struct run*) * sol->width));
	for (int i = 0; i < p->width; i++) {
		sol->col_runs[i] = (struct run*)(malloc(sizeof(struct run) * sol->col_sizes[i]));
	}

	sol->cells = (struct cell**)(malloc(sizeof(struct cell*) * sol->height));
	for (int i = 0; i < p->height; i++) {
		sol->cells[i] = (struct cell*)(malloc(sizeof(struct cell) * sol->width));
		for (int j = 0; j < p->width; j++) {
			sol->cells[i][j].runs = (struct run**)(malloc(sizeof(struct run*) * (sol->row_sizes[i] + sol->col_sizes[j])));
		}
	}

	return sol;
}

void initialize_runs(Puzzle p, Solver sol) {
	// Initialize rows
	for (int i = 0; i < p->height; i++) {
		int size = sol->row_sizes[i];
		sol->row_runs[i][0].s = 0;
		sol->row_runs[i][size - 1].e = p->width - 1;

		//Initialize starting point for each run
		int start = p->row_constraints[i][0].num + 1;
		for (int j = 1; j < size; j++) {
			sol->row_runs[i][j].s = start;
			start += p->row_constraints[i][j].num + 1;
		}

		// Initialize ending point for each run
		int end = p ->width - p->row_constraints[i][size - 1].num - 2;
		for (int j = size - 2; j >= 0; j--) {
			sol->row_runs[i][j].e = end;
			end -= p->row_constraints[i][j].num - 1;
		}

		for (int j = 0; j < size; j++) {
			sol->row_runs[i][j].l = p->row_constraints[i][j].num;
			int runStart = sol->row_runs[i][j].s;
			int runEnd = sol->row_runs[i][j].e + sol->row_runs[i][j].l;
			for (int k = runStart; k < runEnd; k++) {
				// sol->cells[i][k].runs
			}
		}
	}

	// Initialize columns
	for (int i = 0; i < p->width; i++) {
		int size = sol->col_sizes[i];
		sol->col_runs[i][0].s = 0;
		sol->col_runs[i][size - 1].e = p->height - 1;

		// Initialize starting point for each run
		int start = p->col_constraints[i][0].num + 1;
		for (int j = 1; j < size; j++) {
			sol->col_runs[i][j].s = start;
			start += p->col_constraints[i][j].num + 1;
		}

		// Initialize ending point for each run
		int end = p ->height - p->col_constraints[i][size - 1].num - 2;
		for (int j = size - 2; j >= 0; j--) {
			sol->col_runs[i][j].e = end;
			end -= p->col_constraints[i][j].num - 1;
		}

		for (int j = 0; j < size; j++) {
			sol->col_runs[i][j].l = p->col_constraints[i][j].num;
		}
	}

	return;
}

State create_state(Puzzle p, Solution s, Solver sol) {
	int width = p->width;
	int height = p->height;

	State newState = (struct state*)(malloc(sizeof(struct state)));
	Solver solNew = initialize_solver(p);

	solNew->width = width;
	solNew->height = height;
	solNew->row_sizes = sol->row_sizes;
	solNew->col_sizes = sol->col_sizes;

	for (int i = 0; i < height; i++) {
		int size = sol->row_sizes[i];
		for (int j = 0; j < size; j++) {
			solNew->row_runs[i][j] = sol->row_runs[i][j];
		}
	}
	for (int i = 0; i < width; i++) {
		int size = sol->col_sizes[i];
		for (int j = 0; j < size; j++) {
			solNew->col_runs[i][j] = sol->col_runs[i][j];
		}
	}

	newState->sol = solNew;

	Solution sNew = initialize_solution(width, height);
	int size = width * height;
	for (int i = 0; i < size; i++) {
		sNew->data[i] = s->data[i];
	}
}

void solve(Puzzle p, Solution s) {
	s->mark_unknown();

	Solver sol = initialize_solver(p);

	initialize_runs(p, sol);

	int progress = true;
	int iterations = 0;
	while (progress) {
		progress = false;
		printf("Iteration %d\n", iterations);

		for (int i = 0; i < p->height; i++) {
			int size = sol->row_sizes[i];

			// Rule 1.1
			for (int j = 0; j < size; j++) {
				int start = sol->row_runs[i][j].s;
				int end = sol->row_runs[i][j].e;
				int u = end - start + 1 - sol->row_runs[i][j].l;

				for (int k = start + u; k <= end - u; k++) {
					if (s->set(i, k, p->row_constraints[i][j].color)) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = sol->row_runs[i][0].s;
			int lastEnd = sol->row_runs[i][size - 1].e;
			for (int j = 0; j < p->width; j++) {
				if (j < firstStart || j > lastEnd) {
					if (s->set(i, j, EMPTY)) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = sol->row_runs[i][j].e;
				int nextStart = sol->row_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					if (s->set(i, k, EMPTY)) progress = true;
				}
			}

			// Rule 2.1
			for (int j = 1; j < size; j++) {
				int currentStart = sol->row_runs[i][j].s;
				int prevStart = sol->row_runs[i][j-1].s;
				if (currentStart <= prevStart) {
					sol->row_runs[i][j].s = prevStart + sol->row_runs[i][j-1].l + 1;
					progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = sol->row_runs[i][j].e;
				int nextEnd = sol->row_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					sol->row_runs[i][j].e = nextEnd - sol->row_runs[i][j+1].l - 1;
					progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = sol->row_runs[i][j].s;
				int currentEnd = sol->row_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = s->data[i * s->width + currentStart - 1];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						sol->row_runs[i][j].s++;
						progress = true;
					}
				}
				if (currentEnd < s->width - 1) {
					int nextCell = s->data[i * s->width + currentEnd + 1];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						sol->row_runs[i][j].e--;
						progress = true;
					}
				}
			}
		}

		for (int i = 0; i < p->width; i++) {
			int size = sol->col_sizes[i];

			// Rule 1.1
			for (int j = 0; j < size; j++) {
				int start = sol->col_runs[i][j].s;
				int end = sol->col_runs[i][j].e;
				int u = end - start + 1 - sol->col_runs[i][j].l;

				for (int k = start + u; k <= end - u; k++) {
					if (s->set(k, i, p->col_constraints[i][j].color)) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = sol->col_runs[i][0].s;
			int lastEnd = sol->col_runs[i][size - 1].e;
			for (int j = 0; j < p->height; j++) {
				if (j < firstStart || j > lastEnd) {
					if (s->set(j, i, EMPTY)) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = sol->col_runs[i][j].e;
				int nextStart = sol->col_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					if (s->set(k, i, EMPTY)) progress = true;
				}
			}

			// Rule 2.1
			for (int j = 1; j < size; j++) {
				int currentStart = sol->col_runs[i][j].s;
				int prevStart = sol->col_runs[i][j-1].s;
				if (currentStart <= prevStart) {
					sol->col_runs[i][j].s = prevStart + sol->col_runs[i][j-1].l + 1;
					progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = sol->col_runs[i][j].e;
				int nextEnd = sol->col_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					sol->col_runs[i][j].e = nextEnd - sol->col_runs[i][j+1].l - 1;
					progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = sol->col_runs[i][j].s;
				int currentEnd = sol->col_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = s->data[(currentStart - 1) * s->width + i];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						sol->col_runs[i][j].s++;
						progress = true;
					}
				}
				if (currentEnd < s->width - 1) {
					int nextCell = s->data[(currentEnd + 1) * s->width + i];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						sol->col_runs[i][j].e--;
						progress = true;
					}
				}
			}
		}

		s->print_solution();

		iterations++;
	}

	// s->fill_unknown();
	return;
}
