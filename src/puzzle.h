#ifndef __PUZZLE_H__
#define __PUZZLE_H__

struct constraint {
	int num;
	int color;
};

struct puzzle {
	int width;
	int height;

	constraint** row_constraints;
	int* row_sizes;
	constraint** col_constraints;
	int* col_sizes;

	puzzle(int w, int h, int* r_sizes, int* c_sizes,
		int* r_nums, int* r_colors, int* c_nums, int* c_colors) {
		width = w;
		height = h;
		row_sizes = r_sizes;
		col_sizes = c_sizes;
		row_constraints = new constraint*[height];
		col_constraints = new constraint*[width];

		int counter = 0;
		for (int row = 0; row < height; row++) {
			int size = r_sizes[row];
			row_constraints[row] = new constraint[size];

			for (int i = 0; i < size; i++) {
				row_constraints[row][i].num = r_nums[counter];
				row_constraints[row][i].color = r_colors[counter];
				counter++;
			}
		}

		counter = 0;
		for (int col = 0; col < width; col++) {
			int size = c_sizes[col];
			col_constraints[col] = new constraint[size];

			for (int i = 0; i < size; i++) {
				col_constraints[col][i].num = c_nums[counter];
				col_constraints[col][i].color = c_colors[counter];
				counter++;
			}
		}
	}
};

using Puzzle = puzzle*;

struct solution {
	int width;
	int height;
	int* data;

	solution(int w, int h) {
		width = w;
		height = h;
		data = new int[width * height];
		for (int i = 0; i < width * height; i++) {
			data[i] = 0;
		}
	}

	void mark_unknown() {
		for (int i = 0; i < width * height; i++) {
			data[i] = -1;
		}
		return;
	}
};

using Solution = solution*;

bool check_solution(Puzzle p, Solution s);

#endif
