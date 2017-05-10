#include "solver.h"
#include <omp.h>
#include "../util/CycleTimer.h"
#include <stdio.h>
#include <stdlib.h>
#include <stack>

#define NUM_THREADS 8

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
    solv->solved_rows = (bool*)(calloc(solv->height, sizeof(bool)));

    solv->col_runs = (struct run**)(malloc(sizeof(struct run*) * solv->width));
    for (int i = 0; i < p->width; i++) {
        solv->col_runs[i] = (struct run*)(malloc(sizeof(struct run) * solv->col_sizes[i]));
    }
    solv->solved_cols = (bool*)(calloc(solv->width, sizeof(bool)));

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
    //  free(solv->row_runs[i]);
    // }
    // for (int i = 0; i < solv->width; i++) {
    //  free(solv->col_runs[i]);
    // }

    free(solv->row_runs);
    free(solv->col_runs);
    free(solv->row_sizes);
    free(solv->col_sizes);
    free(solv);

    return;
}

int global_height; //declare here cus im fucking lazy 
double pragmaT; 
State create_state(Puzzle p, Solution solu, Solver solv) {
    int width = p->width;
    int height = p->height;

    State newState = (struct state*)(malloc(sizeof(struct state)));
    Solver solvNew = initialize_solver(p);

    solvNew->width = width;
    solvNew->height = height;
    global_height = height; 
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

int local_set(int data[], int row, int color) {
        if (row < 0 || row >= global_height) 
            return OUTOFBOUNDS;

        if (data[row] == UNKNOWN) {
            data[row] = color;
            return PROGRESS;
        }
        else if (data[row] == color) {
            return SAME;
        }
        else {
            /* printf("********************CONFLICT********************\n"); */
            return CONFLICT;
        }
}

bool row_solved(int row, Solution s, Puzzle p, Solver solv) {
    if (solv->solved_rows[row]) return true;

    int col = 0;
    int numConstraints = p->row_sizes[row];
    for (int i = 0; i < numConstraints; i++) {
        
        // Ignore leading empty spaces
        while (col < s->width && s->data[row * s->width + col] == EMPTY) {
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

    // Check trailing empty spaces
    while (col < s->width) {
        if (s->data[row * s->width + col] != 0)
            return false;
        col++;
    }

    solv->solved_rows[row] = true;

    return true;
}

bool col_solved(int col, Solution s, Puzzle p, Solver solv) {
    if (solv->solved_cols[col]) return true;

    int row = 0;
    int numConstraints = p->col_sizes[col];
    for (int i = 0; i < numConstraints; i++) {
        
        // Ignore leading empty spaces
        while (row < s->height && s->data[row * s->width + col] == EMPTY) {
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

    // Check trailing empty spaces
    while (row < s->height) {
        if (s->data[row * s->width + col] != 0)
            return false;
        row++;
    }

    solv->solved_cols[col] = true;

    return true;
}

bool solve(Puzzle p, Solution sol) {
    sol->mark_unknown();

    Solver solvStart = initialize_solver(p);

    initialize_runs(p, solvStart);
    
    State stStart = create_state(p, sol, solvStart);
    stStart->row = 0;
    stStart->run_index = -1;

    std::stack<State> stateStack;
    stateStack.push(stStart);

    bool solved = false;

    int width = p->width;
    int height = p->height;

    while (!stateStack.empty()) {
        State st = stateStack.top();
        stateStack.pop();

        Solver solv = st->solv;
        Solution solu = st->solu;

        bool progress = true;
        bool conflict = false;
        int iterations = 0;

        double start = CycleTimer::currentSeconds();
        #pragma omp parallel num_threads(NUM_THREADS)
        {
            int tid = omp_get_thread_num();

            while (progress && !conflict) {
                // printf("Iteration %d\n", iterations);

                #pragma omp barrier
                #pragma single
                {
                    progress = false;
                }
                /* #pragma omp barrier */
                // ======================
                // ======== ROWS ========
                // ======================
                for (int i = tid; i < height; i += NUM_THREADS) {
                    if (!row_solved(i, solu, p, solv)) {
                /* #pragma omp for */  
                /* for (int i = 0; i < height; i++){ */
                    int size = solv->row_sizes[i];
                    int local[width]; //get local copy
                    bool lconflict = false;
                    bool lprogress = false; 
                    for(int j = 0; j < width; j++){
                        local[j] = solu->data[i*width + j];
                    }

                    // ---- PART 1 ----
                    // Rule 1.1
                    for (int j = 0; j < size; j++) {
                        int start = solv->row_runs[i][j].s;
                        int end = solv->row_runs[i][j].e;
                        int u = end - start + 1 - solv->row_runs[i][j].l;

                        for (int k = start + u; k <= end - u; k++) {
                            int status = local_set(local, k, p->row_constraints[i][j].color);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 1.2
                    int firstStart = solv->row_runs[i][0].s;
                    int lastEnd = solv->row_runs[i][size - 1].e;
                    for (int j = 0; j < width; j++) {
                        if (j < firstStart || j > lastEnd) {
                            int status = local_set(local, j, EMPTY);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }
                    for (int j = 0; j < size - 1; j++) {
                        int currentEnd = solv->row_runs[i][j].e;
                        int nextStart = solv->row_runs[i][j+1].s;
                        for (int k = currentEnd + 1; k < nextStart; k++) {
                            int status = local_set(local, k, EMPTY);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 1.3
                    for (int j = 0; j < size; j++) {
                        int start = solv->row_runs[i][j].s;
                        int end = solv->row_runs[i][j].e;
                        int color = p->row_constraints[i][j].color;

                        if (start >= 1 && local[start] == color) {
                            bool length1 = true;
                            for (int k = 0; k < j; k++) {
                                if (solv->row_runs[i][k].s <= start && start <= solv->row_runs[i][k].e && solv->row_runs[i][k].l != 1)
                                    length1 = false;
                            }

                            if (length1) {
                                int status = local_set(local, start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }

                        if (end <= width - 2 && local[end] == color) {
                            bool length1 = true;
                            for (int k = j + 1; k < size; k++) {
                                if (solv->row_runs[i][k].s <= end && end <= solv->row_runs[i][k].e && solv->row_runs[i][k].l != 1)
                                    length1 = false;
                            }

                            if (length1) {
                                int status = local_set(local, end + 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    // Rule 1.4
                    for (int j = 1; j < width - 1; j++) {
                        if (local[(j - 1)] > 0 && local[j] == UNKNOWN && local[(j + 1)] > 0) {
                            int len = 1;
                            for (int k = j - 1; k >= 0 && local[k] > 0; k--) { len++; }
                            for (int k = j + 1; k < width && local[k] > 0; k++) { len++; }

                            int maxLen = 0;
                            for (int k = 0; k < size; k++) {
                                if (solv->row_runs[i][k].s <= j - 1 && solv->row_runs[i][k].e >= j + 1 && solv->row_runs[i][k].l > maxLen)
                                    maxLen = solv->row_runs[i][k].l;
                            }

                            if (len > maxLen) {
                                int status = local_set(local, j, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    // Rule 1.5
                    for (int j = 1; j < width; j++) {
                        int color = 1;

                        if ((local[(j - 1)] == EMPTY || local[(j - 1)] == UNKNOWN) && local[j] > 0) {
                            int minLen = width + 1;
                            for (int k = 0; k < size; k++) {
                                if (solv->row_runs[i][k].s <= j && solv->row_runs[i][k].e >= j && solv->row_runs[i][k].l < minLen)
                                    minLen = solv->row_runs[i][k].l;
                            }

                            if (minLen <= width) {
                                int empty = j - 1;
                                for (; empty >= j - minLen && empty >= 0 && local[empty] != EMPTY; empty--) {}
                                if (empty >= j - minLen + 1) {
                                    for (int k = j + 1; k <= empty + minLen; k++) {
                                        int status = local_set(local, k, color);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                }

                                empty = j + 1;
                                for (; empty <= j + minLen && empty < width && local[empty] != EMPTY; empty++) {}
                                if (empty <= j + minLen - 1) {
                                    for (int k = empty - minLen; k <= j - 1; k++) {
                                        int status = local_set(local, k, color);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                }
                            }


                            int len = -1;
                            int start = j;
                            int end = j;
                            for (; start >= 0 && local[start] > 0; start--) { len++; }
                            start++;
                            for (; end < width && local[end] > 0; end++) { len++; }
                            end--;

                            bool sameLen = true;
                            for (int k = 0; k < size; k++) {
                                if (solv->row_runs[i][k].s <= j && solv->row_runs[i][k].e >= j && solv->row_runs[i][k].l != len)
                                    sameLen = false;
                            }

                            if (sameLen) {
                                int status = local_set(local, start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                                status = local_set(local, end + 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    // ---- PART 2 ----
                    // Rule 2.1
                    for (int j = 1; j < size; j++) {
                        int currentStart = solv->row_runs[i][j].s;
                        int prevStart = solv->row_runs[i][j-1].s;
                        if (currentStart <= prevStart) {
                            int status = solv->row_runs[i][j].start(prevStart + solv->row_runs[i][j-1].l + 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }
                    for (int j = 0; j < size - 1; j++) {
                        int currentEnd = solv->row_runs[i][j].e;
                        int nextEnd = solv->row_runs[i][j+1].e;
                        if (currentEnd >= nextEnd) {
                            int status = solv->row_runs[i][j].end(nextEnd - solv->row_runs[i][j+1].l - 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 2.2
                    for (int j = 0; j < size; j++) {
                        int currentStart = solv->row_runs[i][j].s;
                        int currentEnd = solv->row_runs[i][j].e;
                        if (currentStart > 0) {
                            int prevCell = local[currentStart - 1];
                            if (prevCell != EMPTY && prevCell != UNKNOWN) {
                                int status = solv->row_runs[i][j].start(solv->row_runs[i][j].s + 1);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                        if (currentEnd < solu->width - 1) {
                            int nextCell = local[currentEnd + 1];
                            if (nextCell != EMPTY && nextCell != UNKNOWN) {
                                int status = solv->row_runs[i][j].end(solv->row_runs[i][j].e - 1);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
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
                            if (local[k] == color) {
                                segEnd = k;
                            }
                            else {
                                if (segEnd - segStart + 1 > len) {
                                    if (segEnd <= solv->row_runs[i][j-1].e && segStart < solv->row_runs[i][j+1].s) {
                                        int status = solv->row_runs[i][j].start(segEnd + 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                    if (segStart >= solv->row_runs[i][j+1].s && segEnd > solv->row_runs[i][j-1].e) {
                                        int status = solv->row_runs[i][j].end(segStart - 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
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
                        // Rule 3.1
                        int prevEnd = j == 0 ? -1 : solv->row_runs[i][j-1].e;
                        int nextStart = j == size - 1 ? width : solv->row_runs[i][j+1].s;
                        int startCell = prevEnd + 1;
                        for (; startCell < nextStart && local[startCell] <= 0; startCell++) {}
                        int endCell = nextStart - 1;
                        for (; endCell > prevEnd && local[endCell] <= 0; endCell--) {}

                        int u = solv->row_runs[i][j].l - (endCell - startCell + 1);
                        if (startCell <= endCell && u >= 0) {
                            for (int k = startCell + 1; k < endCell; k++) {
                                int status = local_set(local, k, p->row_constraints[i][j].color);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }

                            if (startCell - u > solv->row_runs[i][j].s) {
                                int status = solv->row_runs[i][j].start(startCell - u);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (endCell + u < solv->row_runs[i][j].e) {
                                int status = solv->row_runs[i][j].end(endCell + u);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                        //Rule 3.2
                        int start = solv->row_runs[i][j].s;
                        int end = solv->row_runs[i][j].e;
                        int len = solv->row_runs[i][j].l;
                        int segLen = 0;
                        int index = start;
                        for (int k = start; k <= end; k++) {
                            if (local[k] != EMPTY) {
                                segLen++;
                            }
                            if (local[k] == EMPTY || k == end) {
                                if (segLen >= len) {
                                    int status = solv->row_runs[i][j].start(index);
                                    if (status == CONFLICT) lconflict = true;
                                    if (status == PROGRESS) lprogress = true;
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
                            if (local[k] != EMPTY) {
                                segLen++;
                            }
                            if (local[k] == EMPTY || k == start) {
                                if (segLen >= len) {
                                    int status = solv->row_runs[i][j].end(index);
                                    if (status == CONFLICT) lconflict = true;
                                    if (status == PROGRESS) lprogress = true;
                                }
                                else {
                                    segLen = 0;
                                    index = k - 1;
                                }
                            }
                        }
                    

                        //Rule 3.3.1 
                        int color = p->row_constraints[i][j].color;
                        if (local[start] == color &&
                            (j == 0 || solv->row_runs[i][j-1].e < start)) {
                            for (int k = start + 1; k <= start + len - 1; k++) {
                                int status = local_set(local, k, color);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (start - 1 >= 0) {
                                int status = local_set(local, start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (start + len <= width - 1) {
                                int status = local_set(local, start + len, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            solv->row_runs[i][j].e = start + len - 1;
                            if (j < size - 1 && solv->row_runs[i][j+1].s <= solv->row_runs[i][j].e) {
                                int status = solv->row_runs[i][j+1].start(solv->row_runs[i][j].e + 2);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (j > 0 && solv->row_runs[i][j-1].e >= start - 1){
                                int status = solv->row_runs[i][j-1].end(start - 2);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }

                        // Rule 3.3-2
                        int black = start;
                        for (; black < end && local[black] != color; black++) {}
                        int empty = black;
                        for (; empty <= end && local[empty] != EMPTY; empty++) {}
                        if ((j == 0 || start > solv->row_runs[i][j-1].e) &&
                            empty < end && empty > black) {
                            int status = solv->row_runs[i][j].end(empty - 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }


                        if (j == 0 || start > solv->row_runs[i][j-1].e) {
                            int black = start;
                            for (; black < end && local[black] != color; black++) {}

                            int index = black;
                            for (; index <= end && local[index] == color; index++) {}
                            
                            index++;
                            for (int k = index; k <= end; k++) {
                                if (local[k] != color || k == end) {
                                    if ((k - 1) - black + 1 > len) {
                                        int status = solv->row_runs[i][j].end(index - 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                        k = end + 1;
                                    }
                                    index = k + 1;
                                }
                            }
                        }
                    }

                    // Write to global memory
                    for(int j = 0; j < width; j++){
                        solu->set(i,j,local[j]);  
                    }
                    if(lconflict) conflict = true;
                    if(lprogress) progress = true;

                    }
                }

                #pragma omp barrier

                // =========================
                // ======== COLUMNS ========
                // =========================
                if (!conflict)
                {
                for (int i = tid; i < width; i += NUM_THREADS) {
                    if (!col_solved(i, solu, p, solv)) {
                //#pragma omp for  
                /* for (int i = 0; i < width; i++) { */
                    int local_col[height]; //get local copy
                    bool lconflict = false;
                    bool lprogress = false; 
                    for(int j = 0; j < height; j++){
                        local_col[j] = solu->data[j*width + i];
                    }
                    int size = solv->col_sizes[i];

                    // ---- PART 1 ----
                    // Rule 1.1
                    for (int j = 0; j < size; j++) {
                        int start = solv->col_runs[i][j].s;
                        int end = solv->col_runs[i][j].e;
                        int u = end - start + 1 - solv->col_runs[i][j].l;

                        for (int k = start + u; k <= end - u; k++) {
                            int status = local_set(local_col,k,p->col_constraints[i][j].color);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 1.2
                    int firstStart = solv->col_runs[i][0].s;
                    int lastEnd = solv->col_runs[i][size - 1].e;
                    for (int j = 0; j < height; j++) {
                        if (j < firstStart || j > lastEnd) {
                            int status = local_set(local_col,j, EMPTY);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }
                    for (int j = 0; j < size - 1; j++) {
                        int currentEnd = solv->col_runs[i][j].e;
                        int nextStart = solv->col_runs[i][j+1].s;
                        for (int k = currentEnd + 1; k < nextStart; k++) {
                            int status = local_set(local_col,k, EMPTY);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 1.3
                    for (int j = 0; j < size; j++) {
                        int start = solv->col_runs[i][j].s;
                        int end = solv->col_runs[i][j].e;
                        int color = p->col_constraints[i][j].color;

                        if (start >= 1 && local_col[start] == color) {
                            bool length1 = true;
                            for (int k = 0; k < j; k++) {
                                if (solv->col_runs[i][k].s <= start && start <= solv->col_runs[i][k].e && solv->col_runs[i][k].l != 1)
                                    length1 = false;
                            }

                            if (length1) {
                                int status = local_set(local_col,start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }

                        if (end <= height - 2 && local_col[end] == color) {
                            bool length1 = true;
                            for (int k = j + 1; k < size; k++) {
                                if (solv->col_runs[i][k].s <= end && end <= solv->col_runs[i][k].e && solv->col_runs[i][k].l != 1)
                                    length1 = false;
                            }

                            if (length1) {
                                int status = local_set(local_col,end + 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    // Rule 1.4
                    for (int j = 1; j < height - 1; j++) {
                        if (local_col[(j - 1)] > 0 && local_col[j] == UNKNOWN && local_col[(j + 1)] > 0) {
                            int len = 1;
                            for (int k = j - 1; k >= 0 && local_col[k] > 0; k--) { len++; }
                            for (int k = j + 1; k < height && local_col[k] > 0; k++) { len++; }

                            int maxLen = 0;
                            for (int k = 0; k < size; k++) {
                                if (solv->col_runs[i][k].s <= j - 1 && solv->col_runs[i][k].e >= j + 1 && solv->col_runs[i][k].l > maxLen)
                                    maxLen = solv->col_runs[i][k].l;
                            }

                            if (len > maxLen) {
                                int status = local_set(local_col,j, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    // Rule 1.5
                    for (int j = 1; j < height; j++) {
                        int color = 1;

                        if ((local_col[(j - 1)] == EMPTY || local_col[(j - 1) ] == UNKNOWN) && local_col[j] > 0) {
                            int minLen = height + 1;
                            for (int k = 0; k < size; k++) {
                                if (solv->col_runs[i][k].s <= j && solv->col_runs[i][k].e >= j && solv->col_runs[i][k].l < minLen)
                                    minLen = solv->col_runs[i][k].l;
                            }

                            if (minLen <= height) {
                                int empty = j - 1;
                                for (; empty >= j - minLen && empty >= 0 && local_col[empty] != EMPTY; empty--) {}
                                if (empty >= j - minLen + 1) {
                                    for (int k = j + 1; k <= empty + minLen; k++) {
                                        int status = local_set(local_col,k, color);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                }

                                empty = j + 1;
                                for (; empty <= j + minLen && empty < height && local_col[empty] != EMPTY; empty++) {}
                                if (empty <= j + minLen - 1) {
                                    for (int k = empty - minLen; k <= j - 1; k++) {
                                        int status = local_set(local_col,k, color);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                }
                            }


                            int len = -1;
                            int start = j;
                            int end = j;
                            for (; start >= 0 && local_col[start ] > 0; start--) { len++; }
                            start++;
                            for (; end < height && local_col[end] > 0; end++) { len++; }
                            end--;

                            bool sameLen = true;
                            for (int k = 0; k < size; k++) {
                                if (solv->col_runs[i][k].s <= j && solv->col_runs[i][k].e >= j && solv->col_runs[i][k].l != len)
                                    sameLen = false;
                            }

                            if (sameLen) {
                                int status = local_set(local_col,start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                                status = local_set(local_col,end + 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    }

                    
                    // ---- PART 2 ----
                    // Rule 2.1
                    for (int j = 1; j < size; j++) {
                        int currentStart = solv->col_runs[i][j].s;
                        int prevStart = solv->col_runs[i][j-1].s;
                        if (currentStart <= prevStart) {
                            int status = solv->col_runs[i][j].start(prevStart + solv->col_runs[i][j-1].l + 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }
                    for (int j = 0; j < size - 1; j++) {
                        int currentEnd = solv->col_runs[i][j].e;
                        int nextEnd = solv->col_runs[i][j+1].e;
                        if (currentEnd >= nextEnd) {
                            int status = solv->col_runs[i][j].end(nextEnd - solv->col_runs[i][j+1].l - 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    }

                    // Rule 2.2
                    for (int j = 0; j < size; j++) {
                        int currentStart = solv->col_runs[i][j].s;
                        int currentEnd = solv->col_runs[i][j].e;
                        if (currentStart > 0) {
                            int prevCell = local_col[(currentStart - 1)];
                            if (prevCell != EMPTY && prevCell != UNKNOWN) {
                                int status = solv->col_runs[i][j].start(solv->col_runs[i][j].s + 1);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                        if (currentEnd < solu->height - 1) {
                            int nextCell = local_col[(currentEnd + 1)];
                            if (nextCell != EMPTY && nextCell != UNKNOWN) {
                                int status = solv->col_runs[i][j].end(solv->col_runs[i][j].e - 1);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
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
                            if (local_col[k] == color) {
                                segEnd = k;
                            }
                            else {
                                if (segEnd - segStart + 1 > len) {
                                    if (segEnd <= solv->col_runs[i][j-1].e && segStart < solv->col_runs[i][j+1].s) {
                                        int status = solv->col_runs[i][j].start(segEnd + 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                    }
                                    if (segStart >= solv->col_runs[i][j+1].s && segEnd > solv->col_runs[i][j-1].e) {
                                        int status = solv->col_runs[i][j].end(segStart - 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
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
                        for (; startCell < nextStart && local_col[startCell] <= 0; startCell++) {}
                        int endCell = nextStart - 1;
                        for (; endCell > prevEnd && local_col[endCell] <= 0; endCell--) {}

                        int u = solv->col_runs[i][j].l - (endCell - startCell + 1);
                        if (startCell <= endCell && u >= 0) {
                            for (int k = startCell + 1; k < endCell; k++) {
                                int status = local_set(local_col,k, p->col_constraints[i][j].color);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }

                            if (startCell - u > solv->col_runs[i][j].s) {
                                int status = solv->col_runs[i][j].start(startCell - u);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (endCell + u < solv->col_runs[i][j].e) {
                                int status = solv->col_runs[i][j].end(endCell + u);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }
                    // Rule 3.2 
                        int start = solv->col_runs[i][j].s;
                        int end = solv->col_runs[i][j].e;
                        int len = solv->col_runs[i][j].l;
                        int segLen = 0;
                        int index = start;
                        for (int k = start; k <= end; k++) {
                            if (local_col[k ] != EMPTY) {
                                segLen++;
                            }
                            if (local_col[k] == EMPTY || k == end) {
                                if (segLen >= len) {
                                    int status = solv->col_runs[i][j].start(index);
                                    if (status == CONFLICT) lconflict = true;
                                    if (status == PROGRESS) lprogress = true;
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
                            if (local_col[k ] != EMPTY) {
                                segLen++;
                            }
                            if (local_col[k ] == EMPTY || k == start) {
                                if (segLen >= len) {
                                    int status = solv->col_runs[i][j].end(index);
                                    if (status == CONFLICT) lconflict = true;
                                    if (status == PROGRESS) lprogress = true;
                                }
                                else {
                                    segLen = 0;
                                    index = k - 1;
                                }
                            }
                        }
                    //Rule 3.3-1 
                        int color = p->col_constraints[i][j].color;
                        if (local_col[start] == color &&
                            (j == 0 || solv->col_runs[i][j-1].e < start)) {
                            for (int k = start + 1; k <= start + len - 1; k++) {
                                int status = local_set(local_col,k, color);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (start - 1 >= 0) {
                                int status = local_set(local_col,start - 1, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (start + len <= height - 1) {
                                int status = local_set(local_col,start + len, EMPTY);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            solv->col_runs[i][j].e = start + len - 1;
                            if (j < size - 1 && solv->col_runs[i][j+1].s <= solv->col_runs[i][j].e) {
                                int status = solv->col_runs[i][j+1].start(solv->col_runs[i][j].e + 2);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                            if (j > 0 && solv->col_runs[i][j-1].e >= start - 1){
                                int status = solv->col_runs[i][j-1].end(start - 2);
                                if (status == CONFLICT) lconflict = true;
                                if (status == PROGRESS) lprogress = true;
                            }
                        }

                    // Rule 3.3-2
                        int black = start;
                        for (; black < end && local_col[black] != color; black++) {}
                        int empty = black;
                        for (; empty <= end && local_col[empty ] != EMPTY; empty++) {}
                        if ((j == 0 || start > solv->col_runs[i][j-1].e) &&
                            empty < end && empty > black) {
                            int status = solv->col_runs[i][j].end(empty - 1);
                            if (status == CONFLICT) lconflict = true;
                            if (status == PROGRESS) lprogress = true;
                        }
                    // Rule 3.3-3
                        if (j == 0 || start > solv->col_runs[i][j-1].e) {
                            int black = start;
                            for (; black < end && local_col[black] != color; black++) {}

                            int index = black;
                            for (; index <= end && local_col[index] == color; index++) {}
                            
                            index++;
                            for (int k = index; k <= end; k++) {
                                if (local_col[k] != color || k == end) {
                                    if ((k - 1) - black + 1 > len) {
                                        int status = solv->col_runs[i][j].end(index - 2);
                                        if (status == CONFLICT) lconflict = true;
                                        if (status == PROGRESS) lprogress = true;
                                        k = end + 1;
                                    }
                                    index = k + 1;
                                }
                            }
                        }
                    }

                    // Write to global memory
                    for(int j = 0; j < height; j++){
                        solu->set(j,i,local_col[j]);  
                    }
                    if(lconflict) conflict = true;
                    if(lprogress) progress = true;
                    }
                }
                }
                #pragma omp barrier

                // solu->print_solution();
                
                // #pragma omp single
                // {
                //  iterations++;
                // }
            }
        }

        double ref = CycleTimer::currentSeconds() - start;
        pragmaT += ref;

        if (conflict) {
            // printf("Conflict found\n");
            free_state(st);
        }
        else if (filled(solu)) {
            solved = true;

            for (int i = 0; i < sol->width * sol->height; i++) {
                sol->data[i] = st->solu->data[i];
            }
            free_state(st);

            while (!stateStack.empty()) {
                free_state(stateStack.top());
                stateStack.pop();
            }
        }   // Do DFS if not solved
        else {
             // printf("========================Starting DFS========================\n"); 

            int row = st->row;
            int runIndex = st->run_index + 1;
            while (row_solved(row, solu, p, solv)) {
                row++;
                runIndex = 0;
            }
            if (runIndex >= p->row_sizes[row]) {
                row++;
                runIndex = 0;
            }
            if (row >= height) {
                // Failed to solve on this branch
            }
            else {
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

                    if (!conflict) {
                        stateStack.push(newSt);
                    }
                    else {
                        free_state(newSt);
                    }
                }
            }
        }

    }

    printf("TOTAL PRAGMA TIME: %.7f\n", pragmaT);

    return solved;
}

bool filled(Solution solu) {
    int size = solu->width * solu->height;
    for (int i = 0; i < size; i++) {
        if (solu->data[i] == UNKNOWN)
            return false;
    }

    return true;
}
