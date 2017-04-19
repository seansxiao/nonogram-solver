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

	solv->cells = (struct cell**)(malloc(sizeof(struct cell*) * solv->height));
	for (int i = 0; i < p->height; i++) {
		solv->cells[i] = (struct cell*)(malloc(sizeof(struct cell) * solv->width));
		for (int j = 0; j < p->width; j++) {
			solv->cells[i][j].runs = (struct run**)(malloc(sizeof(struct run*) * (solv->row_sizes[i] + solv->col_sizes[j])));
		}
	}

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

State create_state(Puzzle p, Solution solu, Solver solv) {
	int width = p->width;
	int height = p->height;

	State newState = (struct state*)(malloc(sizeof(struct state)));
	Solver solNew = initialize_solver(p);

	solNew->width = width;
	solNew->height = height;
	solNew->row_sizes = solv->row_sizes;
	solNew->col_sizes = solv->col_sizes;

	for (int i = 0; i < height; i++) {
		int size = solv->row_sizes[i];
		for (int j = 0; j < size; j++) {
			solNew->row_runs[i][j] = solv->row_runs[i][j];
		}
	}
	for (int i = 0; i < width; i++) {
		int size = solv->col_sizes[i];
		for (int j = 0; j < size; j++) {
			solNew->col_runs[i][j] = solv->col_runs[i][j];
		}
	}

	newState->solv = solNew;

	Solution sNew = initialize_solution(width, height);
	int size = width * height;
	for (int i = 0; i < size; i++) {
		sNew->data[i] = solu->data[i];
	}

	return newState;
}

void solve(Puzzle p, Solution solu) {
	solu->mark_unknown();

	Solver solv = initialize_solver(p);

	initialize_runs(p, solv);

	int progress = true;
	int iterations = 0;
	while (progress) {
		progress = false;
		printf("Iteration %d\n", iterations);

		for (int i = 0; i < p->height; i++) {
			int size = solv->row_sizes[i];

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
			for (int j = 0; j < p->width; j++) {
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
            /* // Rule 1.3 */
            /* for(int j = 0; j < size; j++){ */
            /*     //start case */
        		/* int start = solv->row_runs[i][j].s; */
            /*     if((start-1) >= 0 && solu->data[i * solu->width + start] > 0){ */
            /*         bool len1 = true; */
            /*         for(int k = 0; k < j; k++){ */
            /*             int diffstart = solv->row_runs[i][k].s; */
            /*             int diffend = solv->row_runs[i][k].e; */
            /*             if(solv->row_runs[i][k].l != 1 || !(diffstart <= start && diffend >= start)){ */
            /*                 len1 = false; */
            /*                 break; */
            /*             } */
            /*         } */
            /*         if(len1){if(solu->set(i,start-1,EMPTY)) progress = true;} */ 
            /*     } */
            /*     //end case */ 
				/* int end = solv->row_runs[i][j].e; */
            /*     if((end+1) < p->width && solu->data[i * solu->width + end] > 0){ */
            /*         bool len1 = true; */
            /*         for(int k = j+1; k < size; k++){ */
            /*             int diffstart = solv->row_runs[i][k].s; */
            /*             int diffend = solv->row_runs[i][k].e; */
            /*             if(solv->row_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)){ */
            /*                 len1 = false; */
            /*                 break; */
            /*             } */
            /*         } */
            /*         if(len1){if(solu->set(i,end+1,EMPTY)) progress = true;} */ 
            /*     } */
            /* } */
            
            /* //Rule 1.4 */
            /* int start_start = -1, start_end = -1;//Ends are not inclusive */ 
            /* int end_start = -1, end_end = -1;//Ends are not inclusive */
            /* int lower_run = 0; //lower bound of runs we have to check for overlap */ 
            /* for(int j = 0; j < p->width; j++){ */
            /*     int cell = solu->data[i*solu->width + j]; */
            /*     if(cell > 0){ */
            /*         if(start_start == -1) start_start = j; */
            /*         if(start_end != -1) end_start = j; */ 
            /*     } */
            /*     else{ */
            /*         if(cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j; */
            /*         else if(end_start != -1 && end_end == -1){//Found black segment - unknown - black segment */
            /*             end_end = j; */ 
            /*             int startlen = start_end - start_start; */
            /*             int endlen = end_end - end_start; */ 
            /*             int totallen = startlen + 1 + endlen; */ 
            /*             int targetcell = start_end; */ 
            /*             //Find any runs that overlap targetcell */
            /*             int max = -1; */
            /*             int min = size; */
            /*             for(int k = lower_run; k < size; k++){ */
            /*                 int runstart = solv->row_runs[i][k].s; */
            /*                 int runend = solv->row_runs[i][k].e; */
            /*                 int runlen = solv->row_runs[i][k].l; */
            /*                 if(runstart <= targetcell && runend >= targetcell){ */
            /*                     max = std::max(max,runlen); */
            /*                     min = std::min(min,k); */ 
            /*                 } */
            /*                 else if (min < size) break; */ 
            /*             } */
            /*             if(max < totallen){ if(solu->set(i,targetcell,EMPTY)) progress = true;} */
            /*             lower_run = min; */ 
            /*             if(cell == UNKNOWN){start_start = end_start, start_end = end_end, end_start = -1, end_end = -1;}//start = end */ 
            /*         } */
            /*         else{ */
            /*             start_start = -1, start_end = -1, end_start = -1, end_end = -1; */
            /*         } */
            /*     } */
            /* } */
            /* // Rule 1.5 */
            int prevEmpty = -1;  
            for(int j = 1; j < p -> width; j++){
                int index = i*solu->width + j;
                int lower_run = 0; 
                if(solu->data[index-1] == EMPTY) prevEmpty = j-1; 
                if(solu->data[index] > 0 && solu->data[index-1] <= 0){

                    int minlen = p -> width;
                    int minindex = size; 
                    for(int k = lower_run; k < size; k++){
                        int runstart = solv->row_runs[i][k].s;
                        int runend = solv->row_runs[i][k].e;
                        int runlen = solv->row_runs[i][k].l;
                        if(runstart <= j && runend >= j){
                            minindex = std::min(minindex,k); 
                            minlen = std::min(minlen, runlen); 
                        }
                        else if (minlen < p->width) break;
                        lower_run = minindex; 
                    }
                    if(prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)){
                        //Color each cell in between
                        for(int k = prevEmpty; k < j; k++){
                            if(solu->set(i,k,1)) progress = true; 
                        }
                    }
                    /* //Find prevAfter */
                    /* int prevAfter = -1; */
                    /* for(int k = j+1; k <= j+minlen-1; k++){ */
                    /*     if(solu->data[i*p->width + k] == EMPTY){prevAfter = k; break;} */
                    /* } */
                    /* if(prevAfter != -1){ */
                    /*      for(int k = j+1; k < prevAfter; k++){ */
                    /*         if(solu->set(i,k,1)) progress = true; */ 
                    /*     } */
                    /* } */
                }
            }
            

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

			// Rule 3.1
			for (int j = 1; j < size - 1; j++) {
				int prevEnd = solv->row_runs[i][j-1].e;
				int nextStart = solv->row_runs[i][j+1].s;
				int startCell = prevEnd + 1;
				for (; startCell < nextStart && solu->data[i * solu->width + startCell] == UNKNOWN; startCell++) {}
				int endCell = nextStart - 1;
				for (; endCell > prevEnd && solu->data[i * solu->width + endCell] == UNKNOWN; endCell--) {}
				
				for (int k = startCell; k < endCell; k++) {
					if (solu->set(i, k, p->row_constraints[i][j].color)) progress = true;
				}

				int u = solv->row_runs[i][j].l - (endCell - startCell + 1);
				if (startCell < endCell && u > 0) {
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
		}

//-------------------- Column -------------------------

		for (int i = 0; i < p->width; i++) {
			int size = solv->col_sizes[i];

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
			for (int j = 0; j < p->height; j++) {
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
            /* // Rule 1.3 */
            /* for(int j = 0; j < size; j++){ */
            /*     //start case */
        		/* int start = solv->col_runs[i][j].s; */
            /*     if((start-1) >= 0 && solu->data[solu->width*start + i] > 0){ */
            /*         bool len1 = true; */
            /*         for(int k = 0; k < j; k++){ */
            /*             int diffstart = solv->col_runs[i][k].s; */
            /*             int diffend = solv->col_runs[i][k].e; */
            /*             if(solv->col_runs[i][k].l != 1 || !(diffstart <= start && diffend >= start)){ */
            /*                 len1 = false; */
            /*                 break; */
            /*             } */
            /*         } */
            /*         if(len1){if(solu->set(start-1,i,EMPTY)) progress = true;} */
            /*     } */
            /*     //end case */ 
				/* int end = solv->col_runs[i][j].e; */
            /*     if((end+1) < p->height && solu->data[solu->width*end + i] > 0){ */
            /*         bool len1 = true; */
            /*         for(int k = j+1; k < size; k++){ */
            /*             int diffstart = solv->col_runs[i][k].s; */
            /*             int diffend = solv->col_runs[i][k].e; */
            /*             if(solv->col_runs[i][k].l != 1 || !(diffstart <= end && diffend >= end)){ */
            /*                 len1 = false; */
            /*                 break; */
            /*             } */
            /*         } */
            /*         if(len1){if(solu->set(end+1,i,EMPTY)) progress = true;} */
            /*     } */
            /* } */

            /* //Rule 1.4 */
            /* int start_start = -1, start_end = -1;//Ends are not inclusive */ 
            /* int end_start = -1, end_end = -1;//Ends are not inclusive */
            /* int lower_run = 0; //lower bound of runs we have to check for overlap */ 
            /* for(int j = 0; j < p->height; j++){ */
            /*     int cell = solu->data[j*solu->width + i]; */
            /*     if(cell > 0){ */
            /*         if(start_start == -1) start_start = j; */
            /*         else if(start_end != -1) end_start = j; */ 
            /*     } */
            /*     else{ */
            /*         if(cell == UNKNOWN && start_start != -1 && start_end == -1) start_end = j; */
            /*         if(end_start != -1 && end_end == -1){//Found black segment - unknown - black segment */
            /*             end_end = j; */ 
            /*             int startlen = start_end - start_start; */
            /*             int endlen = end_end - end_start; */ 
            /*             int totallen = startlen + 1 + endlen; */ 
            /*             int targetcell = start_end; */ 
            /*             //Find any runs that overlap targetcell */
            /*             int max = -1; */
            /*             int min = size; */
            /*             for(int k = lower_run; k < size; k++){ */
            /*                 int runstart = solv->col_runs[i][k].s; */
            /*                 int runend = solv->col_runs[i][k].e; */
            /*                 int runlen = solv->col_runs[i][k].l; */
            /*                 if(runstart <= targetcell && runend >= targetcell){ */
            /*                     max = std::max(max,runlen); */
            /*                     min = std::min(min,k); */ 
            /*                 } */
            /*                 else if(min < size) break;//If we can find one and then fail */  
            /*             } */
            /*             if(max < totallen){ if(solu->set(targetcell,i,EMPTY)) progress = true;} */
            /*             lower_run = min; */ 
            /*             if(cell == UNKNOWN){start_start = end_start, start_end = end_end, end_start = -1, end_end = -1;}//start = end */ 
            /*         } */ 
            /*          else{ */
            /*             start_start = -1, start_end = -1, end_start = -1, end_end = -1; */
            /*         } */
            /*     } */
            /* } */

            /* // Rule 1.5 */
            /* int prevEmpty = -1; */  
            /* for(int j = 1; j < p -> height; j++){ */
            /*     int prevIndex = (j-1)*solu->width + i; */
            /*     int index = j*solu->width + i; */
            /*     int lower_run = 0; */ 
            /*     if(solu->data[prevIndex] == EMPTY) prevEmpty = j-1; */ 
            /*     if(solu->data[index] > 0 && solu->data[prevIndex] <= 0){ */

            /*         int minlen = p -> height; */
            /*         int minindex = size; */ 
            /*         for(int k = lower_run; k < size; k++){ */
            /*             int runstart = solv->col_runs[i][k].s; */
            /*             int runend = solv->col_runs[i][k].e; */
            /*             int runlen = solv->col_runs[i][k].l; */
            /*             if(runstart <= j && runend >= j){ */
            /*                 minindex = std::min(minindex,k); */ 
            /*                 minlen = std::min(minlen, runlen); */ 
            /*             } */
            /*             else if (minlen < p->height) break; */
            /*             lower_run = minindex; */ 
            /*         } */
            /*         if(prevEmpty != -1 && prevEmpty >= (j-minlen+1) && prevEmpty <= (j-1)){ */
            /*             //Color each cell in between */
            /*             for(int k = prevEmpty; k < j; k++){ */
            /*                 if(solu->set(k,i,1)) progress = true; */ 
            /*             } */
            /*         } */
            /*         //Find prevAfter */
            /*         int prevAfter = -1; */
            /*         for(int k = j+1; k <= j+minlen-1; k++){ */
            /*             if(solu->data[k*p->width + i] == EMPTY){prevAfter = k; break;} */
            /*         } */
            /*         if(prevAfter != -1){ */
            /*              for(int k = j+1; k < prevAfter; k++){ */
            /*                 if(solu->set(k,i,1)) progress = true; */ 
            /*             } */
            /*         } */       
            /*     } */
            /* } */
            
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
		}

		solu->print_solution();
        printf("start: %d, end: %d", solv->row_runs[0][0].s, solv->row_runs[0][0].e);

		iterations++;
	}

	// Do DFS if not solved
	if (!solved(solu)) {
		printf("Starting DFS\n");
		State newSt = create_state(p, solu, solv);
		for (int i = 0; i < solv->height; i++) {
			int size = solv->row_sizes[i];
			for (int j = 0; j < size; j++) {

			}
		}
	}

	// solu->fill_unknown();
	return;
}

bool solved(Solution solu) {
	int size = solu->width * solu->height;
	for (int i = 0; i < size; i++) {
		if (solu->data[i] == UNKNOWN)
			return false;
	}

	return true;
}
