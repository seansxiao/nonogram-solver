#include "solver.h"

#include <stdio.h>
#include <stdlib.h>

void solve(Puzzle p, Solution s) {
	s->mark_unknown();

	for (int row = 0; row < p->height; row++) {
		int numConstraints = p->row_sizes[row];
		int* runStarts = new int[numConstraints];
		int* runEnds = new int[numConstraints];
		runStarts[0] = 0;
		runEnds[numConstraints-1] = p->width - 1;


		// Set begin and end points for each run
		int startCounter = p->row_constraints[row][0].num + 1;
		for (int i = 1; i < numConstraints; i++) {
			runStarts[i] = startCounter;
			startCounter += p->row_constraints[row][i].num + 1;
		}

		int endCounter = p->width - p->row_constraints[row][numConstraints-1].num - 2;
		for (int i = numConstraints - 2; i >= 0; i--) {
			runEnds[i] = endCounter;
			endCounter -= p->row_constraints[row][i].num - 1;
		}

		for (int i = 0; i < numConstraints; i++) {
			int u = runEnds[i] - runStarts[i] + 1 - p->row_constraints[row][i].num;

			for (int j = runStarts[i] + u; j <= runEnds[i] - u; j++) {
				s->data[row * s->width + j] = 1;
			}
		}
	}

	for (int col = 0; col < p->width; col++) {
		int numConstraints = p->col_sizes[col];
		int* runStarts = new int[numConstraints];
		int* runEnds = new int[numConstraints];
		runStarts[0] = 0;
		runEnds[numConstraints-1] = p->height - 1;


		// Set begin and end points for each run
		int startCounter = p->col_constraints[col][0].num + 1;
		for (int i = 1; i < numConstraints; i++) {
			runStarts[i] = startCounter;
			startCounter += p->col_constraints[col][i].num + 1;
		}

		int endCounter = p->height - p->col_constraints[col][numConstraints-1].num - 2;
		for (int i = numConstraints - 2; i >= 0; i--) {
			runEnds[i] = endCounter;
			endCounter -= p->col_constraints[col][i].num - 1;
		}

		for (int i = 0; i < numConstraints; i++) {
			int u = runEnds[i] - runStarts[i] + 1 - p->col_constraints[col][i].num;

			for (int j = runStarts[i] + u; j <= runEnds[i] - u; j++) {
				s->data[j * s->width + col] = 1;
			}
		}
	}

	return;
}