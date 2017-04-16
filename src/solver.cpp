#include "solver.h"

#include <stdio.h>
#include <stdlib.h>

bool check_solution(puzzle* p, solution* s) {
	if (p->width != s->width || p->height != s->height)
		return false;

	for (int row = 0; row < p->height; row++) {
		
		int numConstraints = p->row_sizes[row];
		for (int i = 0; i < numConstraints; i++) {
			
			// Ignore leading empty spaces
			int col = 0;
			while (col < s->width && s->data[row * s->width + col] == -1) {
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

	return true;
}