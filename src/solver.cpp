#include "solver.h"

#include <stdio.h>
#include <stdlib.h>

bool check_solution(puzzle* p, solution* s) {
	if (p->width != s->width || p->height != s->height)
		return false;

	// Check row constraints
	for (int row = 0; row < p->height; row++) {
		
		int col = 0;
		int numConstraints = p->row_sizes[row];
		for (int i = 0; i < numConstraints; i++) {
			
			// Ignore leading empty spaces
			while (col < s->width && s->data[row * s->width + col] == 0) {
				col++;
			}

			int color = p->row_constraints[row][i].color;
			int num = p->row_constraints[row][i].num;
			for (int j = 0; j < num; j++) {
				// Out of bounds
				if (col >= s->width)
					return false;
				if (s->data[row * s->width + col] != color)
					return false;
				col++;
			}
		}
	}

	// Check column constraints
	for (int col = 0; col < p->width; col++) {
		
		int row = 0;
		int numConstraints = p->col_sizes[col];
		for (int i = 0; i < numConstraints; i++) {
			
			// Ignore leading empty spaces
			while (row < s->height && s->data[row * s->width + col] == 0) {
				row++;
			}

			int color = p->col_constraints[col][i].color;
			int num = p->col_constraints[col][i].num;
			for (int j = 0; j < num; j++) {
				// Out of bounds
				if (row >= s->height)
					return false;
				if (s->data[row * s->width + col] != color)
					return false;
				row++;
			}
		}
	}

	return true;
}