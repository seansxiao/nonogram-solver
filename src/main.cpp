#include <stdio.h>
#include <stdlib.h>

#include "puzzle.h"
#include "solver.h"

int main() {
	// int w = 8;
	// int h = 11;
	// int r_sizes[11] = {1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1};
	// int c_sizes[8] = {1, 1, 1, 2, 2, 1, 1, 1};
	// int r_nums[13] = {0, 4, 6, 2, 2, 2, 2, 6, 4, 2, 2, 2, 0};
	// int r_colors[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	// int c_nums[10] = {0, 9, 9, 2, 2, 2, 2, 4, 4, 0};
	// int c_colors[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	// int w = 4;
	// int h = 3;
	// int r_sizes[3] = {1, 2, 1};
	// int c_sizes[4] = {1, 2, 1, 1};
	// int r_nums[4] = {2, 1, 1, 4};
	// int r_colors[4] = {1, 1, 1, 1};
	// int c_nums[5] = {2, 1, 1, 3, 1};
	// int c_colors[5] = {1, 1, 1, 1, 1};

	int w = 15;
	int h = 15;
	int r_sizes[15] = {1, 2, 2, 4, 3, 4, 3, 3, 4, 2, 2, 3, 3, 3, 5};
	int c_sizes[15] = {2, 2, 2, 3, 2, 4, 3, 3, 2, 4, 5, 2, 3, 3, 2};
	int r_nums[44] = {3, 2, 5, 4, 6, 2, 3, 2, 1, 1, 6, 1, 2, 3, 4, 2, 5, 1, 1, 6, 1, 4, 3, 4, 1, 1, 1, 5, 5, 1, 2, 3, 1, 3, 3, 1, 2, 2, 1, 2, 3, 3, 1, 1};
	int r_colors[44] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	int c_nums[42] = {3, 1, 6, 1, 2, 3, 2, 6, 3, 7, 5, 4, 2, 3, 1, 2, 1, 1, 2, 2, 3, 5, 7, 4, 2, 1, 4, 3, 1, 1, 1, 1, 5, 3, 2, 2, 6, 1, 1, 2, 1, 1};
	int c_colors[42] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

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

	Solution s = initialize_solution(w, h);
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
	solve(&p, s);

	bool correct = check_solution(&p, s);
	if (correct)
		printf("CORRECTNESS PASSED\n");
	else
		printf("CORRECTNESS FAILED\n");

	return 0;
}
