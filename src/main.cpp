#include <stdio.h>
#include <stdlib.h>

#include "solver.h"

int main() {
	int w = 8;
	int h = 11;
	int r_sizes[11] = {1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1};
	int c_sizes[8] = {1, 1, 1, 2, 2, 1, 1, 1};
	int r_nums[13] = {0, 4, 6, 2, 2, 2, 2, 6, 4, 2, 2, 2, 0};
	int r_colors[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	int c_nums[10] = {0, 9, 9, 2, 2, 2, 2, 4, 4, 0};
	int c_colors[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	puzzle p = puzzle(w, h, r_sizes, c_sizes,
		r_nums, r_colors, c_nums, c_colors);

	// Print column constraints
	for (int i = 0; i < p.width; i++) {
		for (int j = 0; j < p.col_sizes[i]; j++) {
			printf("%d ", p.col_constraints[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	// Print row constraints
	for (int i = 0; i < p.height; i++) {
		for (int j = 0; j < p.row_sizes[i]; j++) {
			printf("%d ", p.row_constraints[i][j]);
		}
		printf("\n");
	}

	return 0;
}
