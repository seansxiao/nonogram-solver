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
	bool solved = solve_helper(p, st);

	for (int i = 0; i < solu->width * solu->height; i++) {
		solu->data[i] = st->solu->data[i];
	}

	free_state(st);

	printf("Solution:\n");
	solu->print_solution();
	if (solved)
		printf("Solved!\n");
	else
		printf("UNABLE TO SOLVE\n");

	return solved;
}

bool solve_helper(Puzzle p, State st) {
	int width = p->width;
	int height = p->height;

	Solver solv = st->solv;
	Solution solu = st->solu;

	int progress = true;
	bool conflict = false;
	int iterations = 0;
	// #pragma omp parallel
	// {
	while (progress && !conflict) {
		progress = false;
		// printf("Iteration %d\n", iterations);

		// ======================
		// ======== ROWS ========
		// ======================
		#pragma omp parallel for
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
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->row_runs[i][0].s;
			int lastEnd = solv->row_runs[i][size - 1].e;
			for (int j = 0; j < width; j++) {
				if (j < firstStart || j > lastEnd) {
					int status = solu->set(i, j, EMPTY);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->row_runs[i][j].e;
				int nextStart = solv->row_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					int status = solu->set(i, k, EMPTY);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.3
			for (int j = 0; j < size; j++) {
				int start = solv->row_runs[i][j].s;
				int end = solv->row_runs[i][j].e;
				int color = p->row_constraints[i][j].color;

				if (start >= 1 && solu->data[i * width + start] == color) {
					bool length1 = true;
					for (int k = 0; k < j; k++) {
						if (solv->row_runs[i][k].s <= start && start <= solv->row_runs[i][k].e && solv->row_runs[i][k].l != 1)
							length1 = false;
					}

					if (length1) {
						int status = solu->set(i, start - 1, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}

				if (end <= width - 2 && solu->data[i * width + end] == color) {
					bool length1 = true;
					for (int k = j + 1; k < size; k++) {
						if (solv->row_runs[i][k].s <= end && end <= solv->row_runs[i][k].e && solv->row_runs[i][k].l != 1)
							length1 = false;
					}

					if (length1) {
						int status = solu->set(i, end + 1, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.4
			for (int j = 1; j < width - 1; j++) {
				if (solu->data[i * width + (j - 1)] > 0 && solu->data[i * width + j] == UNKNOWN && solu->data[i * width + (j + 1)] > 0) {
					int len = 1;
					for (int k = j - 1; k >= 0 && solu->data[i * width + k] > 0; k--) { len++; }
					for (int k = j + 1; k < width && solu->data[i * width + k] > 0; k++) { len++; }

					int maxLen = 0;
					for (int k = 0; k < size; k++) {
						if (solv->row_runs[i][k].s <= j - 1 && solv->row_runs[i][k].e >= j + 1 && solv->row_runs[i][k].l > maxLen)
							maxLen = solv->row_runs[i][k].l;
					}

					if (len > maxLen) {
						int status = solu->set(i, j, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.5
			for (int j = 1; j < width; j++) {
				int color = 1;

				if ((solu->data[i * width + (j - 1)] == EMPTY || solu->data[i * width + (j - 1)] == UNKNOWN) && solu->data[i * width + j] > 0) {
					int minLen = width + 1;
					for (int k = 0; k < size; k++) {
						if (solv->row_runs[i][k].s <= j && solv->row_runs[i][k].e >= j && solv->row_runs[i][k].l < minLen)
							minLen = solv->row_runs[i][k].l;
					}

					if (minLen <= width) {
						int empty = j - 1;
						for (; empty >= j - minLen && empty >= 0 && solu->data[i * width + empty] != EMPTY; empty--) {}
						if (empty >= j - minLen + 1) {
							for (int k = j + 1; k <= empty + minLen; k++) {
								int status = solu->set(i, k, color);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
						}

						empty = j + 1;
						for (; empty <= j + minLen && empty < width && solu->data[i * width + empty] != EMPTY; empty++) {}
						if (empty <= j + minLen - 1) {
							for (int k = empty - minLen; k <= j - 1; k++) {
								int status = solu->set(i, k, color);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
						}
					}


					int len = -1;
					int start = j;
					int end = j;
					for (; start >= 0 && solu->data[i * width + start] > 0; start--) { len++; }
					start++;
					for (; end < width && solu->data[i * width + end] > 0; end++) { len++; }
					end--;

					bool sameLen = true;
					for (int k = 0; k < size; k++) {
						if (solv->row_runs[i][k].s <= j && solv->row_runs[i][k].e >= j && solv->row_runs[i][k].l != len)
							sameLen = false;
					}

					if (sameLen) {
						int status = solu->set(i, start - 1, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
						status = solu->set(i, end + 1, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.3
			/*for (int j = 0; j < size; j++) {
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
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
				// End case 
				int end = solv->row_runs[i][j].e;
                if((end+1) < p->width && solu->data[i * solu->width + end] > 0){
                    bool len1 = true;
                    for(int k = j+1; k < size; k++){
                        int diffstart = solv->row_runs[i][k].s;
                        int diffend = solv->row_runs[i][k].e;
                        if(solv->row_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)){
                            len1 = false;
                            break;
                        }
                    }
                    if(len1){
                        int status = solu->set(i,end+1,EMPTY);
        				if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
                    }
                }
            }
            
            // Rule 1.4
            int start_start = -1, start_end = -1;//Ends are not inclusive 
            int end_start = -1, end_end = -1;//Ends are not inclusive
            int lower_run = 0; //lower bound of runs we have to check for overlap 
            for(int j = 0; j < p->width; j++){
                int cell = solu->data[i*solu->width + j];
                if(cell > 0){
                    if(start_start == -1) start_start = j;
                    if(start_end != -1) end_start = j; 
                }
                else{
                    if(cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j;
                    else if(end_start != -1 && end_end == -1){//Found black segment - unknown - black segment
                        end_end = j; 
                        int startlen = start_end - start_start;
                        int endlen = end_end - end_start; 
                        int totallen = startlen + 1 + endlen; 
                        int targetcell = start_end; 
                        //Find any runs that overlap targetcell
                        int max = -1;
                        int min = size;
                        for(int k = 0; k < size; k++){
                            int runstart = solv->row_runs[i][k].s;
                            int runend = solv->row_runs[i][k].e;
                            int runlen = solv->row_runs[i][k].l;
                            if(runstart <= targetcell && runend >= targetcell){
                                max = std::max(max,runlen);
                                min = std::min(min,k); 
                            }
                            //else if (min < size) break; 
                        }
                        if(max < totallen){
                            int status = solu->set(i,targetcell,EMPTY);
                            if(status == CONFLICT) conflict = true;
                            if(status == PROGRESS) progress = true;
                        }
                        lower_run = min; 
                        if(cell == UNKNOWN){start_start = end_start, start_end = end_end, end_start = -1, end_end = -1;}//start = end 
                    }
                    else{
                        start_start = -1, start_end = -1, end_start = -1, end_end = -1;
                    }
                }
            }

           // Rule 1.5
            int prevEmpty = -1;  
            //int lower_run = 0; 
            for(int j = 1; j < p -> width; j++){
                int index = i*solu->width + j;
                if(solu->data[index-1] == EMPTY) prevEmpty = j-1; 
                if(solu->data[index] > 0 && solu->data[index-1] <= 0){

                    int minlen = p -> width;
                    //int minindex = size;
                    for(int k = 0; k < size; k++){
                        int runstart = solv->row_runs[i][k].s;
                        int runend = solv->row_runs[i][k].e;
                        int runlen = solv->row_runs[i][k].l;
                        if(runstart <= j && runend >= j){
                            //minindex = std::min(minindex,k); 
                            minlen = std::min(minlen, runlen); 
                        }
                        //else if (minlen < p->width) break;
                    }
                    //lower_run = minindex; 
                    if(prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)){
                        //Color each cell in between
                        for(int k = j+1; k <= prevEmpty + minlen; k++){
                            int status = solu->set(i,k,1); 
                            if(status == CONFLICT) conflict = true;
                            if(status == PROGRESS) progress = true;
                        }
                    }
                    //Find afterEmpty
                    int afterEmpty = -1;
                    for(int k = j+1; k <= j+minlen-1; k++){
                        if(solu->data[i*p->width + k] == EMPTY){afterEmpty = k; break;}
                    }
                    if(afterEmpty != -1){
                         for(int k = afterEmpty-minlen; k <= j-1; k++){
                            int status = solu->set(i,k,1); 
                            if(status == CONFLICT) conflict = true;
                            if(status == PROGRESS) progress = true;
 
                        }
                    }
                }
            }
            //Rule 1.6
            for(int j = 1; j < p->width; j++){
                int index = i*p->width + j; 
                if(solu->data[index] > 0){
                    int segment = 1;
                    for(int k = j+1; k < p->width; k++){
                        if(solu->data[i*p->width + k] > 0) segment++; 
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
                        int status = solu->set(i,j-1,EMPTY);
                        if(status == CONFLICT) conflict = true;
                        if(status == PROGRESS) progress = true;
 
                        status = solu->set(i,j+segment-1,EMPTY); 
                        // if(status == CONFLICT) conflict = true;
                        if(status == PROGRESS) progress = true;
 
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
					int status = solv->row_runs[i][j].start(prevStart + solv->row_runs[i][j-1].l + 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->row_runs[i][j].e;
				int nextEnd = solv->row_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					int status = solv->row_runs[i][j].end(nextEnd - solv->row_runs[i][j+1].l - 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = solv->row_runs[i][j].s;
				int currentEnd = solv->row_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = solu->data[i * solu->width + currentStart - 1];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						int status = solv->row_runs[i][j].start(solv->row_runs[i][j].s + 1);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
				if (currentEnd < solu->width - 1) {
					int nextCell = solu->data[i * solu->width + currentEnd + 1];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						int status = solv->row_runs[i][j].end(solv->row_runs[i][j].e - 1);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
								int status = solv->row_runs[i][j].start(segEnd + 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
							if (segStart >= solv->row_runs[i][j+1].s && segEnd > solv->row_runs[i][j-1].e) {
								int status = solv->row_runs[i][j].end(segStart - 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
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
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}

					if (startCell - u > solv->row_runs[i][j].s) {
						int status = solv->row_runs[i][j].start(startCell - u);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (endCell + u < solv->row_runs[i][j].e) {
						int status = solv->row_runs[i][j].end(endCell + u);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
						if (segLen >= len) {
							int status = solv->row_runs[i][j].start(index);
							if (status == CONFLICT) conflict = true;
							if (status == PROGRESS) progress = true;
						}
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
						if (segLen >= len) {
							int status = solv->row_runs[i][j].end(index);
							if (status == CONFLICT) conflict = true;
							if (status == PROGRESS) progress = true;
						}
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
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (start - 1 >= 0) {
						int status = solu->set(i, start - 1, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (start + len <= width - 1) {
						int status = solu->set(i, start + len, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					solv->row_runs[i][j].e = start + len - 1;
					if (j < size - 1 && solv->row_runs[i][j+1].s <= solv->row_runs[i][j].e) {
						int status = solv->row_runs[i][j+1].start(solv->row_runs[i][j].e + 2);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (j > 0 && solv->row_runs[i][j-1].e >= start - 1){
						int status = solv->row_runs[i][j-1].end(start - 2);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
					int status = solv->row_runs[i][j].end(empty - 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
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
								int status = solv->row_runs[i][j].end(index - 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
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
		#pragma omp parallel for
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
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.2
			int firstStart = solv->col_runs[i][0].s;
			int lastEnd = solv->col_runs[i][size - 1].e;
			for (int j = 0; j < height; j++) {
				if (j < firstStart || j > lastEnd) {
					int status = solu->set(j, i, EMPTY);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->col_runs[i][j].e;
				int nextStart = solv->col_runs[i][j+1].s;
				for (int k = currentEnd + 1; k < nextStart; k++) {
					int status = solu->set(k, i, EMPTY);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 1.3
			for (int j = 0; j < size; j++) {
				int start = solv->col_runs[i][j].s;
				int end = solv->col_runs[i][j].e;
				int color = p->col_constraints[i][j].color;

				if (start >= 1 && solu->data[start * width + i] == color) {
					bool length1 = true;
					for (int k = 0; k < j; k++) {
						if (solv->col_runs[i][k].s <= start && start <= solv->col_runs[i][k].e && solv->col_runs[i][k].l != 1)
							length1 = false;
					}

					if (length1) {
						int status = solu->set(start - 1, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}

				if (end <= height - 2 && solu->data[end * width + i] == color) {
					bool length1 = true;
					for (int k = j + 1; k < size; k++) {
						if (solv->col_runs[i][k].s <= end && end <= solv->col_runs[i][k].e && solv->col_runs[i][k].l != 1)
							length1 = false;
					}

					if (length1) {
						int status = solu->set(end + 1, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.4
			for (int j = 1; j < height - 1; j++) {
				if (solu->data[(j - 1) * width + i] > 0 && solu->data[j * width + i] == UNKNOWN && solu->data[(j + 1) * width + i] > 0) {
					int len = 1;
					for (int k = j - 1; k >= 0 && solu->data[k * width + i] > 0; k--) { len++; }
					for (int k = j + 1; k < height && solu->data[k * width + i] > 0; k++) { len++; }

					int maxLen = 0;
					for (int k = 0; k < size; k++) {
						if (solv->col_runs[i][k].s <= j - 1 && solv->col_runs[i][k].e >= j + 1 && solv->col_runs[i][k].l > maxLen)
							maxLen = solv->col_runs[i][k].l;
					}

					if (len > maxLen) {
						int status = solu->set(j, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.5
			for (int j = 1; j < height; j++) {
				int color = 1;

				if ((solu->data[(j - 1) * width + i] == EMPTY || solu->data[(j - 1) * width + i] == UNKNOWN) && solu->data[j * width + i] > 0) {
					int minLen = height + 1;
					for (int k = 0; k < size; k++) {
						if (solv->col_runs[i][k].s <= j && solv->col_runs[i][k].e >= j && solv->col_runs[i][k].l < minLen)
							minLen = solv->col_runs[i][k].l;
					}

					if (minLen <= height) {
						int empty = j - 1;
						for (; empty >= j - minLen && empty >= 0 && solu->data[empty * width + i] != EMPTY; empty--) {}
						if (empty >= j - minLen + 1) {
							for (int k = j + 1; k <= empty + minLen; k++) {
								int status = solu->set(k, i, color);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
						}

						empty = j + 1;
						for (; empty <= j + minLen && empty < height && solu->data[empty * width + i] != EMPTY; empty++) {}
						if (empty <= j + minLen - 1) {
							for (int k = empty - minLen; k <= j - 1; k++) {
								int status = solu->set(k, i, color);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
						}
					}


					int len = -1;
					int start = j;
					int end = j;
					for (; start >= 0 && solu->data[start * width + i] > 0; start--) { len++; }
					start++;
					for (; end < height && solu->data[end * width + i] > 0; end++) { len++; }
					end--;

					bool sameLen = true;
					for (int k = 0; k < size; k++) {
						if (solv->col_runs[i][k].s <= j && solv->col_runs[i][k].e >= j && solv->col_runs[i][k].l != len)
							sameLen = false;
					}

					if (sameLen) {
						int status = solu->set(start - 1, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
						status = solu->set(end + 1, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
			}

			// Rule 1.3
			/*for (int j = 0; j < size; j++) {
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
    					if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
				// End case 
				int end = solv->col_runs[i][j].e;
                if((end+1) < p->height && solu->data[solu->width*end + i] > 0){
                    bool len1 = true;
                    for(int k = j+1; k < size; k++){
                        int diffstart = solv->col_runs[i][k].s;
                        int diffend = solv->col_runs[i][k].e;
                        if(solv->col_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)){
                            len1 = false;
                            break;
                        }
                    }
                    if(len1){
                        int status = solu->set(end+1,i,EMPTY);
    					if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
                    }
	
                }
            }

            // Rule 1.4
            int start_start = -1, start_end = -1;//Ends are not inclusive 
            int end_start = -1, end_end = -1;//Ends are not inclusive
            int lower_run = 0; //lower bound of runs we have to check for overlap 
            for(int j = 0; j < p->height; j++){
                int cell = solu->data[j*solu->width + i];
                if(cell > 0){
                    if(start_start == -1) start_start = j;
                    else if(start_end != -1) end_start = j; 
                }
                else{
                    if(cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j;
                    if(end_start != -1 && end_end == -1){//Found black segment - unknown - black segment
                        end_end = j; 
                        int startlen = start_end - start_start;
                        int endlen = end_end - end_start; 
                        int totallen = startlen + 1 + endlen; 
                        int targetcell = start_end; 
                        //Find any runs that overlap targetcell
                        int max = -1;
                        int min = size;
                        for(int k = 0; k < size; k++){
                            int runstart = solv->col_runs[i][k].s;
                            int runend = solv->col_runs[i][k].e;
                            int runlen = solv->col_runs[i][k].l;
                            if(runstart <= targetcell && runend >= targetcell){
                                max = std::max(max,runlen);
                                min = std::min(min,k); 
                            }
                            //else if(min < size) break;//If we can find one and then fail  
                        }
                        if(max < totallen){ 
                            int status = solu->set(targetcell,i,EMPTY);
    				    	if (status == CONFLICT) conflict = true;
						    if (status == PROGRESS) progress = true;
                        }
                        lower_run = min; 
                        if(cell == UNKNOWN){start_start = end_start, start_end = end_end, end_start = -1, end_end = -1;}//start = end 
                    } 
                     else{
                        start_start = -1, start_end = -1, end_start = -1, end_end = -1;
                    }
                }
            }

            //Rule 1.5 
            int prevEmpty = -1;  

            //int lower_run = 0; 
            for(int j = 1; j < p -> height; j++){
                int prevIndex = (j-1)*solu->width + i;
                int index = j*solu->width + i;
                if(solu->data[prevIndex] == EMPTY) prevEmpty = j-1; 
                if(solu->data[index] > 0 && solu->data[prevIndex] <= 0){
                   
                    int minlen = p -> height;
                    int minindex = size; 

                    for(int k = 0; k < size; k++){
                        int runstart = solv->col_runs[i][k].s;
                        int runend = solv->col_runs[i][k].e;
                        int runlen = solv->col_runs[i][k].l;
                        if(runstart <= j && runend >= j){
                            //minindex = std::min(minindex,k); 
                            minlen = std::min(minlen, runlen); 
                        }
                        //else if (minlen < p->height) break;
                        //lower_run = minindex; 
                    }
                    if(prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)){
                        //Color each cell in between
                        for(int k = j+1; k <= prevEmpty+minlen; k++){
                            int status = solu->set(k,i,1);  
    	    				if (status == CONFLICT) conflict = true;
			    			if (status == PROGRESS) progress = true;
	
                        }
                    }
                    //Find afterEmpty
                    int afterEmpty = -1;
                    for(int k = j+1; k <= j+minlen-1; k++){
                        if(solu->data[k*p->width + i] == EMPTY){afterEmpty = k; break;}
                    }
                    if(afterEmpty != -1){
                         for(int k = afterEmpty-minlen; k <= j-1 ; k++){
                            int status = solu->set(k,i,1); 
    				    	if (status == CONFLICT) conflict = true;
						    if (status == PROGRESS) progress = true;
	
                        }
                    }       
                }
            }
            //Rule 1.6
            for(int j = 1; j < p->height; j++){
                int index = j*p->width + i; 
                if(solu->data[index] > 0){
                    int segment = 1;
                    for(int k = j+1; k < p->height; k++){
                        if(solu->data[k*p->width + i] > 0) segment++; 
                        else break; 
                    }   
                    bool samesize = true; 
                    for(int k = 0; k < size; k++){
                        int runstart = solv->col_runs[i][k].s;
                        int runend = solv->col_runs[i][k].e;
                        int runlen = solv->col_runs[i][k].l;
                        if(runstart <= j && runend >= j){
                            samesize = samesize && (runlen == segment); 
                        }
                    }
                    if(samesize){
                        int status = solu->set(j-1,i,EMPTY);
    					if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
	
                        status = solu->set(j+segment-1,i,EMPTY);
    					//if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
	
                        j = j+segment;
                    } 
                }
            }*/


            // ---- PART 2 ----
			// Rule 2.1
			for (int j = 1; j < size; j++) {
				int currentStart = solv->col_runs[i][j].s;
				int prevStart = solv->col_runs[i][j-1].s;
				if (currentStart <= prevStart) {
					int status = solv->col_runs[i][j].start(prevStart + solv->col_runs[i][j-1].l + 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}
			for (int j = 0; j < size - 1; j++) {
				int currentEnd = solv->col_runs[i][j].e;
				int nextEnd = solv->col_runs[i][j+1].e;
				if (currentEnd >= nextEnd) {
					int status = solv->col_runs[i][j].end(nextEnd - solv->col_runs[i][j+1].l - 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
				}
			}

			// Rule 2.2
			for (int j = 0; j < size; j++) {
				int currentStart = solv->col_runs[i][j].s;
				int currentEnd = solv->col_runs[i][j].e;
				if (currentStart > 0) {
					int prevCell = solu->data[(currentStart - 1) * solu->width + i];
					if (prevCell != EMPTY && prevCell != UNKNOWN) {
						int status = solv->col_runs[i][j].start(solv->col_runs[i][j].s + 1);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
				}
				if (currentEnd < solu->height - 1) {
					int nextCell = solu->data[(currentEnd + 1) * solu->width + i];
					if (nextCell != EMPTY && nextCell != UNKNOWN) {
						int status = solv->col_runs[i][j].end(solv->col_runs[i][j].e - 1);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
								int status = solv->col_runs[i][j].start(segEnd + 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
							}
							if (segStart >= solv->col_runs[i][j+1].s && segEnd > solv->col_runs[i][j-1].e) {
								int status = solv->col_runs[i][j].end(segStart - 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
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
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}

					if (startCell - u > solv->col_runs[i][j].s) {
						int status = solv->col_runs[i][j].start(startCell - u);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (endCell + u < solv->col_runs[i][j].e) {
						int status = solv->col_runs[i][j].end(endCell + u);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
						if (segLen >= len) {
							int status = solv->col_runs[i][j].start(index);
							if (status == CONFLICT) conflict = true;
							if (status == PROGRESS) progress = true;
						}
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
						if (segLen >= len) {
							int status = solv->col_runs[i][j].end(index);
							if (status == CONFLICT) conflict = true;
							if (status == PROGRESS) progress = true;
						}
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
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (start - 1 >= 0) {
						int status = solu->set(start - 1, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (start + len <= height - 1) {
						int status = solu->set(start + len, i, EMPTY);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					solv->col_runs[i][j].e = start + len - 1;
					if (j < size - 1 && solv->col_runs[i][j+1].s <= solv->col_runs[i][j].e) {
						int status = solv->col_runs[i][j+1].start(solv->col_runs[i][j].e + 2);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
					}
					if (j > 0 && solv->col_runs[i][j-1].e >= start - 1){
						int status = solv->col_runs[i][j-1].end(start - 2);
						if (status == CONFLICT) conflict = true;
						if (status == PROGRESS) progress = true;
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
					int status = solv->col_runs[i][j].end(empty - 1);
					if (status == CONFLICT) conflict = true;
					if (status == PROGRESS) progress = true;
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
								int status = solv->col_runs[i][j].end(index - 2);
								if (status == CONFLICT) conflict = true;
								if (status == PROGRESS) progress = true;
								k = end + 1;
							}
							index = k + 1;
						}
					}
				}
			}
		}

		// solu->print_solution();

		iterations++;
	}
	// }

	if (conflict) {
		// printf("Conflict found\n");
		return false;
	}

	// Do DFS if not solved
	if (filled(solu)) {
		return true;
	}
	else {
		//printf("========================Starting DFS========================\n");
		State newSt = create_state(p, solu, solv);
		int row = st->row;
		int runIndex = st->run_index + 1;
		if (runIndex >= p->row_sizes[row]) {
			row++;
			runIndex = 0;
		}
		if (row >= height)
			return false;
		// printf("row: %d, run: %d\n", row, runIndex);
		int runStart = st->solv->row_runs[row][runIndex].s;
		int runRight = st->solv->row_runs[row][runIndex].e - st->solv->row_runs[row][runIndex].l + 1;
		for (; runStart <= runRight; runStart++) {
			State newSt = create_state(p, solu, solv);
			newSt->row = row;
			newSt->run_index = runIndex;

			int runEnd = runStart + st->solv->row_runs[row][runIndex].l - 1;

			// Set run at particular location
			newSt->solv->row_runs[row][runIndex].s = runStart;
			newSt->solv->row_runs[row][runIndex].e = runEnd;

			bool conflict = false;
			// Update all cells in the run region
			for (int i = runStart; i <= runEnd; i++) {
				if (newSt->solu->set(row, i, p->row_constraints[row][runIndex].color) == CONFLICT)
					conflict = true;
			}
			if (runStart > 0) {
				if (newSt->solu->set(row, runStart - 1, EMPTY) == CONFLICT)
					conflict = true;
			}
			if (runEnd < width - 1) {
				if (newSt->solu->set(row, runEnd + 1, EMPTY) == CONFLICT)
					conflict = true;
			}

			if (!conflict && solve_helper(p, newSt)) {
				*st = *newSt;
				return true;
			}
			else {
				free_state(newSt);
			}
		}
	}

	return false;
}

bool filled(Solution solu) {
	int size = solu->width * solu->height;
	for (int i = 0; i < size; i++) {
		if (solu->data[i] == UNKNOWN)
			return false;
	}

	return true;
}
