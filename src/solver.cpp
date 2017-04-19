#include "solver.h"

#include <stdio.h>
#include <stdlib.h>

Solver initialize_solver(Puzzle p) {
	solver* solv = (struct solver*)(malloc(sizeof(struct solver)));

	solv->width = p->width;
	solv->height = p->height;
	
	solv->row_sizes = p->row_sizes;
	solv->col_sizes = p->col_sizes;

	solv->row_runs = (struct run**)(malloc(sizeof(struct run*) * solv->height));
	for (int i = 0; i < p->height; i++) {
		solv->row_runs[i] = (struct run*)(malloc(sizeof(struct run) * solv->row_sizes[i]));
	}

	solv->col_runs = (struct run**)(malloc(sizeof(struct run*) * solv->width));
	for (int i = 0; i < p->width; i++) {
		solv->col_runs[i] = (struct run*)(malloc(sizeof(struct run) * solv->col_sizes[i]));
	}

	// solv->cells = (struct cell**)(malloc(sizeof(struct cell*) * solv->height));
	// for (int i = 0; i < p->height; i++) {
	// 	solv->cells[i] = (struct cell*)(malloc(sizeof(struct cell) * solv->width));
	// 	for (int j = 0; j < p->width; j++) {
	// 		solv->cells[i][j].runs = (struct run**)(malloc(sizeof(struct run*) * (solv->row_sizes[i] + solv->col_sizes[j])));
	// 	}
	// }

	return solv;
}

void initialize_runs(Puzzle p, Solver solv) {
	// Initialize rows
	for (int i = 0; i < p->height; i++) {
		int size = solv->row_sizes[i];
		solv->row_runs[i][0].s = 0;
		solv->row_runs[i][size - 1].e = p->width - 1;

		//Initialize starting point for each run
		int start = p->row_constraints[i][0].num + 1;
		for (int j = 1; j < size; j++) {
			solv->row_runs[i][j].s = start;
			start += p->row_constraints[i][j].num + 1;
		}

		// Initialize ending point for each run
		int end = p ->width - p->row_constraints[i][size - 1].num - 2;
		for (int j = size - 2; j >= 0; j--) {
			solv->row_runs[i][j].e = end;
			end -= p->row_constraints[i][j].num - 1;
		}

		for (int j = 0; j < size; j++) {
			solv->row_runs[i][j].l = p->row_constraints[i][j].num;
			int runStart = solv->row_runs[i][j].s;
			int runEnd = solv->row_runs[i][j].e + solv->row_runs[i][j].l;
			for (int k = runStart; k < runEnd; k++) {
				// solv->cells[i][k].runs
			}
		}
	}

	// Initialize columns
	for (int i = 0; i < p->width; i++) {
		int size = solv->col_sizes[i];
		solv->col_runs[i][0].s = 0;
		solv->col_runs[i][size - 1].e = p->height - 1;

		// Initialize starting point for each run
		int start = p->col_constraints[i][0].num + 1;
		for (int j = 1; j < size; j++) {
			solv->col_runs[i][j].s = start;
			start += p->col_constraints[i][j].num + 1;
		}

		// Initialize ending point for each run
		int end = p ->height - p->col_constraints[i][size - 1].num - 2;
		for (int j = size - 2; j >= 0; j--) {
			solv->col_runs[i][j].e = end;
			end -= p->col_constraints[i][j].num - 1;
		}

		for (int j = 0; j < size; j++) {
			solv->col_runs[i][j].l = p->col_constraints[i][j].num;
		}
	}

	return;
}

void free_solver(Solver solv) {
	// for (int i = 0; i < solv->height; i++) {
	// 	free(solv->row_runs[i]);
	// }
	// for (int i = 0; i < solv->width; i++) {
	// 	free(solv->col_runs[i]);
	// }
	free(solv->row_runs);
	free(solv->col_runs);
	free(solv->row_sizes);
	free(solv->col_sizes);
	free(solv);
	return;
}

State create_state(Puzzle p, Solution solu, Solver solv) {
	int width = p->width;
	int height = p->height;

	State newState = (struct state*)(malloc(sizeof(struct state)));
	Solver solvNew = initialize_solver(p);

	solvNew->width = width;
	solvNew->height = height;
	solvNew->row_sizes = solv->row_sizes;
	solvNew->col_sizes = solv->col_sizes;

	for (int i = 0; i < height; i++) {
		int size = solv->row_sizes[i];
		for (int j = 0; j < size; j++) {
			solvNew->row_runs[i][j] = solv->row_runs[i][j];
		}
	}
	for (int i = 0; i < width; i++) {
		int size = solv->col_sizes[i];
		for (int j = 0; j < size; j++) {
			solvNew->col_runs[i][j] = solv->col_runs[i][j];
		}
	}

	newState->solv = solvNew;

	Solution soluNew = initialize_solution(width, height);
	int size = width * height;
	for (int i = 0; i < size; i++) {
		soluNew->data[i] = solu->data[i];
	}

	newState->solu = soluNew;

	return newState;
}

void free_state(State st) {
	// free_solver(st->solv);
	free(st->solu->data);
	free(st->solu);
	free(st);
	return;
}

bool solve(Puzzle p, Solution solu) {
	solu->mark_unknown();

	Solver solv = initialize_solver(p);

	initialize_runs(p, solv);
	
	State st = create_state(p, solu, solv);
	st->row = 0;
	st->run_index = -1;
	solve_helper(p, st);

	for (int i = 0; i < solu->width * solu->height; i++) {
		solu->data[i] = st->solu->data[i];
	}

	free_state(st);

	return false;
}

bool solve_helper(Puzzle p, State st) {
	int width = p->width;
	int height = p->height;

	Solver solv = st->solv;
	Solution solu = st->solu;

	int progress = true;
	int iterations = 0;
	while (progress) {
		progress = false;
		printf("Iteration %d\n", iterations);

		// ======================
		// ======== ROWS ========
		// ======================
		for (int i = 0; i < height; i++) {
			int size = solv->row_sizes[i];

			// ---- PART 1 ----
			// Rule 1.1
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int u = end - start + 1 - solv->row_runs[i][j].l;

				for (int k = start + u; k <= end - u; k++) {
					if (solu->set(i, k, p->row_constraints[i][j].color)) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->row_runs[i][0].s;
			int lastEnd = solv->row_runs[i][size - 1].e;
			for (int j = 0; j < width; j++) {
				if (j < firstStart || j > lastEnd) {
					if (solu->set(i, j, EMPTY)) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->row_runs[i][j].e;
				int nextStart = solv->row_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					if (solu->set(i, k, EMPTY)) progress = true;
				}
			}

			// ---- PART 2 ----
			// Rule 2.1
			for (int j = 1; j < size; j++) {
				int currentStart = solv->row_runs[i][j].s;
				int prevStart = solv->row_runs[i][j-1].s;
				if (currentStart <= prevStart) {
					solv->row_runs[i][j].s = prevStart + solv->row_runs[i][j-1].l + 1;
					progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->row_runs[i][j].e;
				int nextEnd = solv->row_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					solv->row_runs[i][j].e = nextEnd - solv->row_runs[i][j+1].l - 1;
					progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = solv->row_runs[i][j].s;
				int currentEnd = solv->row_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = solu->data[i * solu->width + currentStart - 1];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						solv->row_runs[i][j].s++;
						progress = true;
					}
				}
				if (currentEnd < solu->width - 1) {
					int nextCell = solu->data[i * solu->width + currentEnd + 1];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						solv->row_runs[i][j].e--;
						progress = true;
					}
				}
			}

			// ---- PART 3 ----
			// Rule 3.1
			for (int j = 0; j < size; j++) {
				int prevEnd = j == 0 ? -1 : solv->row_runs[i][j-1].e;
				int nextStart = j == size - 1 ? width : solv->row_runs[i][j+1].s;
				int startCell = prevEnd + 1;
				for (; startCell < nextStart && solu->data[i * solu->width + startCell] <= 0; startCell++) {}
				int endCell = nextStart - 1;
				for (; endCell > prevEnd && solu->data[i * solu->width + endCell] <= 0; endCell--) {}

				int u = solv->row_runs[i][j].l - (endCell - startCell + 1);
				if (startCell < endCell && u >= 0) {
					for (int k = startCell + 1; k < endCell; k++) {
						if (solu->set(i, k, p->row_constraints[i][j].color)) progress = true;
					}

					if (startCell - u > solv->row_runs[i][j].s) {
						solv->row_runs[i][j].s = startCell - u;
						progress = true;
					}
					if (endCell + u < solv->row_runs[i][j].e) {
						solv->row_runs[i][j].e = endCell + u;
						progress = true;
					}
				}
			}

			// Rule 3.2
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int len = solv->row_runs[i][j].l;
				int segLen = 0;
				int index = start;
				for (int k = start; k <= end; k++) {
					if (solu->data[i * solu->width + k] != EMPTY) {
						segLen++;
					}
					else {
						if (segLen >= len)
							solv->row_runs[i][j].s = index;
						else {
							segLen = 0;
							index = k + 1;
						}
					}
				}
				segLen = 0;
				index = end;
				for (int k = end; k >= start; k--) {
					if (solu->data[i * solu->width + k] != EMPTY) {
						segLen++;
					}
					else {
						if (segLen >= len)
							solv->row_runs[i][j].e = index;
						else {
							segLen = 0;
							index = k - 1;
						}
					}
				}
			}
		}

		// =========================
		// ======== COLUMNS ========
		// =========================
		for (int i = 0; i < width; i++) {
			int size = solv->col_sizes[i];

			// ---- PART 1 ----
			// Rule 1.1
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int u = end - start + 1 - solv->col_runs[i][j].l;

				for (int k = start + u; k <= end - u; k++) {
					if (solu->set(k, i, p->col_constraints[i][j].color)) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->col_runs[i][0].s;
			int lastEnd = solv->col_runs[i][size - 1].e;
			for (int j = 0; j < height; j++) {
				if (j < firstStart || j > lastEnd) {
					if (solu->set(j, i, EMPTY)) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->col_runs[i][j].e;
				int nextStart = solv->col_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					if (solu->set(k, i, EMPTY)) progress = true;
				}
			}

			// ---- PART 2 ----
			// Rule 2.1
			for (int j = 1; j < size; j++) {
				int currentStart = solv->col_runs[i][j].s;
				int prevStart = solv->col_runs[i][j-1].s;
				if (currentStart <= prevStart) {
					solv->col_runs[i][j].s = prevStart + solv->col_runs[i][j-1].l + 1;
					progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->col_runs[i][j].e;
				int nextEnd = solv->col_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					solv->col_runs[i][j].e = nextEnd - solv->col_runs[i][j+1].l - 1;
					progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = solv->col_runs[i][j].s;
				int currentEnd = solv->col_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = solu->data[(currentStart - 1) * solu->width + i];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						solv->col_runs[i][j].s++;
						progress = true;
					}
				}
				if (currentEnd < solu->width - 1) {
					int nextCell = solu->data[(currentEnd + 1) * solu->width + i];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						solv->col_runs[i][j].e--;
						progress = true;
					}
				}
			}

			// ---- PART 3 ----
			// Rule 3.1
			for (int j = 0; j < size; j++) {
				int prevEnd = j == 0 ? -1 : solv->col_runs[i][j-1].e;
				int nextStart = j == size - 1 ? height : solv->col_runs[i][j+1].s;
				int startCell = prevEnd + 1;
				for (; startCell < nextStart && solu->data[startCell * solu->width + i] <= 0; startCell++) {}
				int endCell = nextStart - 1;
				for (; endCell > prevEnd && solu->data[endCell * solu->width + i] <= 0; endCell--) {}

				int u = solv->col_runs[i][j].l - (endCell - startCell + 1);
				if (startCell < endCell && u >= 0) {
					for (int k = startCell + 1; k < endCell; k++) {
						if (solu->set(k, i, p->col_constraints[i][j].color)) progress = true;
					}

					if (startCell - u > solv->col_runs[i][j].s) {
						solv->col_runs[i][j].s = startCell - u;
						progress = true;
					}
					if (endCell + u < solv->col_runs[i][j].e) {
						solv->col_runs[i][j].e = endCell + u;
						progress = true;
					}
				}
			}

			// Rule 3.2
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int len = solv->col_runs[i][j].l;
				int segLen = 0;
				int index = start;
				for (int k = start; k <= end; k++) {
					if (solu->data[k * solu->width + i] != EMPTY) {
						segLen++;
					}
					else {
						if (segLen >= len)
							solv->col_runs[i][j].s = index;
						else {
							segLen = 0;
							index = k + 1;
						}
					}
				}
				segLen = 0;
				index = end;
				for (int k = end; k >= start; k--) {
					if (solu->data[k * solu->width + i] != EMPTY) {
						segLen++;
					}
					else {
						if (segLen >= len)
							solv->col_runs[i][j].e = index;
						else {
							segLen = 0;
							index = k - 1;
						}
					}
				}
			}
		}

		solu->print_solution();
		iterations++;
	}

	printf("Solution:\n");
	solu->print_solution();

	// Do DFS if not solved
	if (false && !solved(solu)) {
		printf("Starting DFS\n");
		State newSt = create_state(p, solu, solv);
		newSt->row = st->row;
		newSt->run_index = st->run_index + 1;
		if (newSt->run_index > p->row_sizes[newSt->row]) {
			newSt->row = st->row + 1;
			newSt->run_index = 0;
		}
		int row = newSt->row;
		int runIndex = newSt->run_index;
		if (row >= height)
			return true;
		printf("row: %d, run: %d\n", newSt->row, newSt->run_index);
		int runStart = newSt->solv->row_runs[row][runIndex].s;
		int runRight = newSt->solv->row_runs[row][runIndex].e - newSt->solv->row_runs[row][runIndex].l + 1;
		for (; runStart <= runRight; runStart++) {
			int runEnd = runStart + newSt->solv->row_runs[row][runIndex].l - 1;

			// Set run at particular location
			newSt->solv->row_runs[row][runIndex].s = runStart;
			newSt->solv->row_runs[row][runIndex].e = runEnd;

			// Update all cells in the run region
			for (int i = runStart; i <= runEnd; i++) {
				newSt->solu->set(row, i, p->row_constraints[row][runIndex].color);
			}
			if (runStart > 0)
				newSt->solu->set(row, runStart - 1, EMPTY);
			if (runEnd < width - 1)
				newSt->solu->set(row, runEnd + 1, EMPTY);

			// Verify each column for the run
			for (int i = runStart; i <= runEnd; i++) {

			}

			if (solve_helper(p, newSt))
				return true;
		}
	}
	else {
		return true;
	}

	// solu->fill_unknown();
	return false;
}

bool solved(Solution solu) {
	int size = solu->width * solu->height;
	for (int i = 0; i < size; i++) {
		if (solu->data[i] == UNKNOWN)
			return false;
	}

	return true;
}
