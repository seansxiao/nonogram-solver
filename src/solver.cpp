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
					int status = solu->set(i, k, p->row_constraints[i][j].color);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->row_runs[i][0].s;
			int lastEnd = solv->row_runs[i][size - 1].e;
			for (int j = 0; j < width; j++) {
				if (j < firstStart || j > lastEnd) {
					int status = solu->set(i, j, EMPTY);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->row_runs[i][j].e;
				int nextStart = solv->row_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					int status = solu->set(i, k, EMPTY);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.3
			for (int j = 0; j < size; j++) {
				// Start case
				int start = solv->row_runs[i][j].s;
				if ((start-1) >= 0 && solu->data[i * solu->width + start] > 0) {
					bool len1 = true;
					for (int k = 0; k < j; k++){
						int diffstart = solv->row_runs[i][k].s;
						int diffend = solv->row_runs[i][k].e;
						if (solv->row_runs[i][k].l != 1 || !(diffstart <= start && diffend >= start)) {
							len1 = false;
							break;
						}
					}
					if (len1) {
						int status = solu->set(i, start - 1, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
				}
				// End case 
				int end = solv->row_runs[i][j].e;
				if ((end+1) < p->width && solu->data[i * solu->width + end] > 0) {
					bool len1 = true;
					for (int k = j+1; k < size; k++) {
						int diffstart = solv->row_runs[i][k].s;
						int diffend = solv->row_runs[i][k].e;
						if (solv->row_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)) {
							len1 = false;
							break;
						}
					}
					if(len1) {
						int status = solu->set(i, end + 1, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.4
			int start_start = -1, start_end = -1;	// Ends are not inclusive 
			int end_start = -1, end_end = -1;	// Ends are not inclusive
			int lower_run = 0;	// Lower bound of runs we have to check for overlap 
			for (int j = 0; j < p->width; j++) {
				int cell = solu->data[i*solu->width + j];
				if (cell > 0) {
					if (start_start == -1) start_start = j;
					if (start_end != -1) end_start = j;
				}
				else {
					if (cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j;
					else if (end_start != -1 && end_end == -1){	// Found black segment - unknown - black segment
						end_end = j; 
						int startlen = start_end - start_start;
						int endlen = end_end - end_start;
						int totallen = startlen + 1 + endlen;
						int targetcell = start_end; 
						// Find any runs that overlap target cell
						int max = -1;
						int min = size;
						for (int k = lower_run; k < size; k++) {
							int runstart = solv->row_runs[i][k].s;
							int runend = solv->row_runs[i][k].e;
							int runlen = solv->row_runs[i][k].l;
							if (runstart <= targetcell && runend >= targetcell) {
								max = std::max(max,runlen);
								min = std::min(min,k);
							}
							else if (min < size) break; 
						}
						if (max < totallen) {
							int status = solu->set(i, targetcell, EMPTY);
							if (status == CONFLICT)
								return false;
							if (status == PROGRESS) progress = true;
						}
						lower_run = min;
						if (cell == UNKNOWN) { start_start = end_start, start_end = end_end, end_start = -1, end_end = -1; }	//start = end 
					}
					else {
						start_start = -1, start_end = -1, end_start = -1, end_end = -1;
					}
				}
			}

			// Rule 1.5
			/*int prevEmpty = -1;
			for (int j = 1; j < p -> width; j++) {
				int index = i*solu->width + j;
				int lower_run = 0;
				if (solu->data[index-1] == EMPTY) prevEmpty = j-1;
				if (solu->data[index] > 0 && solu->data[index-1] <= 0) {

					int minlen = p -> width;
					int minindex = size;
					for (int k = lower_run; k < size; k++) {
						int runstart = solv->row_runs[i][k].s;
						int runend = solv->row_runs[i][k].e;
						int runlen = solv->row_runs[i][k].l;
						if (runstart <= j && runend >= j) {
							minindex = std::min(minindex,k);
							minlen = std::min(minlen, runlen);
						}
						// else if (minlen < p->width) break;
					}
					lower_run = minindex; 
					if (prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)) {
						// Color each cell in between
						for (int k = j+1; k <= prevEmpty + minlen; k++) {
							if (solu->set(i,k,1)) progress = true;
						}
					}
					// Find prevAfter
					int prevAfter = -1;
					for (int k = j+1; k <= j+minlen-1; k++) {
						if (solu->data[i*p->width + k] == EMPTY) { prevAfter = k; break; }
					}
					if (prevAfter != -1) {
						for (int k = prevAfter-minlen; k <= j-1; k++) {
							if(solu->set(i,k,1)) progress = true; 
						}
					}
				}
			}*/
            /*//Rule 1.6
            for(int j = 1; j < p->width; j++){
                int index = i*p->width + j; 
                if(solu->data[index] > 0){
                    int segment = 1;
                    for(int k = j+1; k < p->width; k++){
                        if(solu->data[k*p->width + i] > 0) segment++; 
                        else break; 
                    }   
                    bool samesize = true; 
                    for(int k = 0; k < size; k++){
                        int runstart = solv->row_runs[i][k].s;
                        int runend = solv->row_runs[i][k].e;
                        int runlen = solv->row_runs[i][k].l;
                        if(runstart <= j && runend >= j){
                            samesize = samesize && (runlen == segment); 
                        }
                    }
                    if(samesize){
                        if(solu->set(i,j-1,EMPTY)) progress = true;
                        //if(solu->set(i,j+segment,EMPTY)) progress=true;
                        j = j+segment;
                    } 
                }
            }*/

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

			// Rule 2.3
			for (int j = 1; j < size - 1; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int len = solv->row_runs[i][j].l;
				int color = p->row_constraints[i][j].color;
				int segStart = start;
				int segEnd = segStart - 1;
				for (int k = start; k <= end; k++) {
					if (solu->data[i * solu->width + k] == color) {
						segEnd = k;
					}
					else {
						if (segEnd - segStart + 1 > len) {
							if (segEnd <= solv->row_runs[i][j-1].e && segStart < solv->row_runs[i][j+1].s) {
								solv->row_runs[i][j].s = segEnd + 2;
								progress = true;
							}
							if (segStart >= solv->row_runs[i][j+1].s && segEnd > solv->row_runs[i][j-1].e) {
								solv->row_runs[i][j].e = segStart - 2;
								progress = true;
							}
						}
						segStart = k + 1;
						segEnd = segStart - 1;
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
				if (startCell <= endCell && u >= 0) {
					for (int k = startCell + 1; k < endCell; k++) {
						int status = solu->set(i, k, p->row_constraints[i][j].color);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
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
					if (solu->data[i * solu->width + k] == EMPTY || k == end) {
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
					if (solu->data[i * solu->width + k] == EMPTY || k == start) {
						if (segLen >= len)
							solv->row_runs[i][j].e = index;
						else {
							segLen = 0;
							index = k - 1;
						}
					}
				}
			}

			// Rule 3.3-1
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int len = solv->row_runs[i][j].l;
				int color = p->row_constraints[i][j].color;
				if (solu->data[i * solu->width + start] == color &&
					(j == 0 || solv->row_runs[i][j-1].e < start)) {
					for (int k = start + 1; k <= start + len - 1; k++) {
						int status = solu->set(i, k, color);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					if (start - 1 >= 0) {
						int status = solu->set(i, start - 1, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					if (start + len <= width - 1) {
						int status = solu->set(i, start + len, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					solv->row_runs[i][j].e = start + len - 1;
					if (j < size - 1 && solv->row_runs[i][j+1].s <= solv->row_runs[i][j].e) {
						solv->row_runs[i][j+1].s = solv->row_runs[i][j].e + 2;
						progress = true;
					}
					if (j > 0 && solv->row_runs[i][j-1].e >= start - 1){
						solv->row_runs[i][j-1].e = start - 2;
						progress = true;
					}
				}
			}

			// Rule 3.3-2
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int color = p->row_constraints[i][j].color;
				int black = start;
				for (; black < end && solu->data[i * solu->width + black] != color; black++) {}
				int empty = black;
				for (; empty <= end && solu->data[i * solu->width + empty] != EMPTY; empty++) {}
				if ((j == 0 || start > solv->row_runs[i][j-1].e) &&
					empty < end && empty > black) {
					solv->row_runs[i][j].e = empty - 1;
					progress = true;
				}
			}

			// Rule 3.3-3
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int color = p->row_constraints[i][j].color;
				int len = solv->row_runs[i][j].l;

				if (j == 0 || start > solv->row_runs[i][j-1].e) {
					int black = start;
					for (; black < end && solu->data[i * solu->width + black] != color; black++) {}

					int index = black;
					for (; index <= end && solu->data[i * solu->width + index] == color; index++) {}
					
					index++;
					for (int k = index; k <= end; k++) {
						if (solu->data[i * solu->width + k] != color || k == end) {
							if ((k - 1) - black + 1 > len) {
								solv->row_runs[i][j].e = index - 2;
								progress = true;
								k = end + 1;
							}
							index = k + 1;
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
					int status = solu->set(k, i, p->col_constraints[i][j].color);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->col_runs[i][0].s;
			int lastEnd = solv->col_runs[i][size - 1].e;
			for (int j = 0; j < height; j++) {
				if (j < firstStart || j > lastEnd) {
					int status = solu->set(j, i, EMPTY);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->col_runs[i][j].e;
				int nextStart = solv->col_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					int status = solu->set(k, i, EMPTY);
					if (status == CONFLICT)
						return false;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.3
			for (int j = 0; j < size; j++) {
				// Start case
				int start = solv->col_runs[i][j].s;
				if ((start-1) >= 0 && solu->data[solu->width*start + i] > 0) {
					bool len1 = true;
					for (int k = 0; k < j; k++) {
						int diffstart = solv->col_runs[i][k].s;
						int diffend = solv->col_runs[i][k].e;
						if (solv->col_runs[i][k].l != 1 || !(diffstart <= start && diffend >= start)) {
							len1 = false;
							break;
						}
					}
					if (len1) {
						int status = solu->set(start - 1, i, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
				}
				// End case 
				int end = solv->col_runs[i][j].e;
				if ((end+1) < p->height && solu->data[solu->width*end + i] > 0) {
					bool len1 = true;
					for (int k = j+1; k < size; k++) {
						int diffstart = solv->col_runs[i][k].s;
						int diffend = solv->col_runs[i][k].e;
						if (solv->col_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)) {
							len1 = false;
							break;
						}
					}
					if (len1) {
						int status = solu->set(end + 1, i, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.4
			int start_start = -1, start_end = -1;	// Ends are not inclusive 
			int end_start = -1, end_end = -1;	// Ends are not inclusive
			int lower_run = 0;	// Lower bound of runs we have to check for overlap 
			for (int j = 0; j < p->height; j++) {
				int cell = solu->data[j*solu->width + i];
				if (cell > 0) {
					if (start_start == -1) start_start = j;
					else if (start_end != -1) end_start = j; 
				}
				else {
					if (cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j;
					if (end_start != -1 && end_end == -1){	// Found black segment - unknown - black segment
						end_end = j;
						int startlen = start_end - start_start;
						int endlen = end_end - end_start;
						int totallen = startlen + 1 + endlen;
						int targetcell = start_end;
						// Find any runs that overlap targetcell
						int max = -1;
						int min = size;
						for (int k = lower_run; k < size; k++) {
							int runstart = solv->col_runs[i][k].s;
							int runend = solv->col_runs[i][k].e;
							int runlen = solv->col_runs[i][k].l;
							if (runstart <= targetcell && runend >= targetcell) {
								max = std::max(max,runlen);
								min = std::min(min,k);
							}
							else if (min < size) break;	// If we can find one and then fail
						}
						if (max < totallen) {
							int status = solu->set(targetcell, i, EMPTY);
							if (status == CONFLICT)
								return false;
							if (status == PROGRESS) progress = true;
						}
						lower_run = min; 
						if (cell == UNKNOWN) { start_start = end_start, start_end = end_end, end_start = -1, end_end = -1; } // start = end
					}
					else {
						start_start = -1, start_end = -1, end_start = -1, end_end = -1;
					}
				}
			}

			// Rule 1.5 
			/*int prevEmpty = -1;
			for (int j = 1; j < p -> height; j++) {
				int prevIndex = (j-1)*solu->width + i;
				int index = j*solu->width + i;
				int lower_run = 0; 
				if (solu->data[prevIndex] == EMPTY) prevEmpty = j-1;
				if (solu->data[index] > 0 && solu->data[prevIndex] <= 0) {

					int minlen = p -> height;
					int minindex = size; 

					for (int k = lower_run; k < size; k++) {
						int runstart = solv->col_runs[i][k].s;
						int runend = solv->col_runs[i][k].e;
						int runlen = solv->col_runs[i][k].l;
						if (runstart <= j && runend >= j) {
							minindex = std::min(minindex,k);
							minlen = std::min(minlen, runlen);
						}
						// else if (minlen < p->height) break;
						lower_run = minindex;
					}
					if (prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)) {
						// Color each cell in between
						for (int k = j+1; k <= prevEmpty+minlen; k++) {
							if (solu->set(k,i,1)) progress = true;
						}
					}
					// Find prevAfter
					int prevAfter = -1;
					for (int k = j+1; k <= j+minlen-1; k++) {
						if (solu->data[k*p->width + i] == EMPTY) { prevAfter = k; break; }
					}
					if (prevAfter != -1) {
						for (int k = prevAfter-minlen; k <= j-1 ; k++) {
							if (solu->set(k,i,1)) progress = true;
						}
					}
				}
			}*/

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
				if (currentEnd < solu->height - 1) {
					int nextCell = solu->data[(currentEnd + 1) * solu->width + i];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						solv->col_runs[i][j].e--;
						progress = true;
					}
				}
			}

			// Rule 2.3
			for (int j = 1; j < size - 1; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int len = solv->col_runs[i][j].l;
				int color = p->col_constraints[i][j].color;
				int segStart = start;
				int segEnd = segStart - 1;
				for (int k = start; k <= end; k++) {
					if (solu->data[k * solu->width + i] == color) {
						segEnd = k;
					}
					else {
						if (segEnd - segStart + 1 > len) {
							if (segEnd <= solv->col_runs[i][j-1].e && segStart < solv->col_runs[i][j+1].s) {
								solv->col_runs[i][j].s = segEnd + 2;
								progress = true;
							}
							if (segStart >= solv->col_runs[i][j+1].s && segEnd > solv->col_runs[i][j-1].e) {
								solv->col_runs[i][j].e = segStart - 2;
								progress = true;
							}
						}
						segStart = k + 1;
						segEnd = segStart - 1;
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
				if (startCell <= endCell && u >= 0) {
					for (int k = startCell + 1; k < endCell; k++) {
						int status = solu->set(k, i, p->col_constraints[i][j].color);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
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
					if (solu->data[k * solu->width + i] == EMPTY || k == end) {
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
					if (solu->data[k * solu->width + i] == EMPTY || k == start) {
						if (segLen >= len)
							solv->col_runs[i][j].e = index;
						else {
							segLen = 0;
							index = k - 1;
						}
					}
				}
			}

			// Rule 3.3-1
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int len = solv->col_runs[i][j].l;
				int color = p->col_constraints[i][j].color;
				if (solu->data[start * solu->width + i] == color &&
					(j == 0 || solv->col_runs[i][j-1].e < start)) {
					for (int k = start + 1; k <= start + len - 1; k++) {
						int status = solu->set(k, i, color);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					if (start - 1 >= 0) {
						int status = solu->set(start - 1, i, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					if (start + len <= height - 1) {
						int status = solu->set(start + len, i, EMPTY);
						if (status == CONFLICT)
							return false;
						if (status == PROGRESS) progress = true;
					}
					solv->col_runs[i][j].e = start + len - 1;
					if (j < size - 1 && solv->col_runs[i][j+1].s <= solv->col_runs[i][j].e) {
						solv->col_runs[i][j+1].s = solv->col_runs[i][j].e + 2;
						progress = true;
					}
					if (j > 0 && solv->col_runs[i][j-1].e >= start - 1){
						solv->col_runs[i][j-1].e = start - 2;
						progress = true;
					}
				}
			}

			// Rule 3.3-2
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int color = p->col_constraints[i][j].color;
				int black = start;
				for (; black < end && solu->data[black * solu->width + i] != color; black++) {}
				int empty = black;
				for (; empty <= end && solu->data[empty * solu->width + i] != EMPTY; empty++) {}
				if ((j == 0 || start > solv->col_runs[i][j-1].e) &&
					empty < end && empty > black) {
					solv->col_runs[i][j].e = empty - 1;
					progress = true;
				}
			}

			// Rule 3.3-3
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int color = p->col_constraints[i][j].color;
				int len = solv->col_runs[i][j].l;

				if (j == 0 || start > solv->col_runs[i][j-1].e) {
					int black = start;
					for (; black < end && solu->data[black * solu->width + i] != color; black++) {}

					int index = black;
					for (; index <= end && solu->data[index * solu->width + i] == color; index++) {}
					
					index++;
					for (int k = index; k <= end; k++) {
						if (solu->data[k * solu->width + i] != color || k == end) {
							if ((k - 1) - black + 1 > len) {
								solv->col_runs[i][j].e = index - 2;
								progress = true;
								k = end + 1;
							}
							index = k + 1;
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
	if (true || solved(solu)) {
		return true;
	}
	else {
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

			if (solve_helper(p, newSt))
				return true;
		}
	}

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
