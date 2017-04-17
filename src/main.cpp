#include <stdio.h>
#include <stdlib.h>

#include "puzzle.h"
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
	printf("Column constraints:\n");
	for (int i = 0; i < p.width; i++) {
		for (int j = 0; j < p.col_sizes[i]; j++) {
			printf("%d ", p.col_constraints[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("Row constraints:\n");
	// Print row constraints
	for (int i = 0; i < p.height; i++) {
		for (int j = 0; j < p.row_sizes[i]; j++) {
			printf("%d ", p.row_constraints[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	solution s = solution(w, h);
	// for (int i = 0; i < w * h; i++) {
	// 	if (i >= 9 && i <= 12)
	// 		s.data[i] = 1;
	// 	if (i >= 17 && i <= 22)
	// 		s.data[i] = 1;
	// 	if (i >= 25 && i <= 26)
	// 		s.data[i] = 1;
	// 	if (i >= 29 && i <= 30)
	// 		s.data[i] = 1;
	// 	if (i >= 33 && i <= 34)
	// 		s.data[i] = 1;
	// 	if (i >= 37 && i <= 38)
	// 		s.data[i] = 1;
	// 	if (i >= 41 && i <= 46)
	// 		s.data[i] = 1;
	// 	if (i >= 49 && i <= 52)
	// 		s.data[i] = 1;
	// 	if (i >= 57 && i <= 58)
	// 		s.data[i] = 1;
	// 	if (i >= 65 && i <= 66)
	// 		s.data[i] = 1;
	// 	if (i >= 73 && i <= 74)
	// 		s.data[i] = 1;
	// }
	solve(&p, &s);
	// Print solution
	printf("Solution:\n");
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			int color = s.data[i * w + j];
			if (color == 0)
				printf("-");
			else if (color == 1)
				printf("X");
			else if (color == -1)
				printf(".");
			else
				printf(" ");
		}
		printf("\n");
	}
	printf("\n");

	bool correct = check_solution(&p, &s);
	if (correct)
		printf("CORRECTNESS PASSED\n");
	else
		printf("CORRECTNESS FAILED\n");

	return 0;
}
