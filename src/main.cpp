#include <stdio.h>
#include <stdlib.h>
#include <fstream>	// File reading
#include "puzzle.h"
#include "solver.h"
#include <vector>
#include <string.h>
#include <sstream>
#include <iostream>
using namespace std; 

int main(int argc, char *argv[]) {
    if(argc != 2){printf("Expected one file name!\n"); return 0;}

    ifstream myfile;
    string line;
    bool row = true;


    myfile.open(argv[1]);
    int w = 0, h = 0;
    vector<int> r_sizes,r_nums,r_colors,c_sizes, c_nums,c_colors;
    // Export from here: http://webpbn.com/export.cgi 
    // Option: Format for Andrew Makhorin's pbnsol nonogram solver. Not for multi-color puzzles.
    while(getline(myfile,line))
    {
        if (line[0] == '*') { continue; }
        if (line[0] == '&') { row = false; continue; }
        if (row) {
            h++;
            string buf; // Have a buffer string
            stringstream ss(line); // Insert the string into a stream
            vector<string> tokens; // Create vector to hold our words
            while (ss >> buf)
                tokens.push_back(buf);

            r_sizes.push_back(tokens.size());
            for(size_t i = 0; i < tokens.size(); i++) {
                r_nums.push_back(stoi(tokens[i])); 
                r_colors.push_back(1);
            }

        }
        else {
            w++;
            string buf;	// Have a buffer string
            stringstream ss(line);	// Insert the string into a stream
            vector<string> tokens;	// Create vector to hold our words
            while (ss >> buf)
                tokens.push_back(buf);

            c_sizes.push_back(tokens.size());
            for (size_t i = 0; i < tokens.size(); i++) {
                c_nums.push_back(stoi(tokens[i])); 
                c_colors.push_back(1);
            }

        }
    }
    myfile.close(); 

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

	/* int w = 15; */
	/* int h = 15; */
	/* int r_sizes[15] = {1, 2, 2, 4, 3, 4, 3, 3, 4, 2, 2, 3, 3, 3, 5}; */
	/* int c_sizes[15] = {2, 2, 2, 3, 2, 4, 3, 3, 2, 4, 5, 2, 3, 3, 2}; */
	/* int r_nums[44] = {3, 2, 5, 4, 6, 2, 3, 2, 1, 1, 6, 1, 2, 3, 4, 2, 5, 1, 1, 6, 1, 4, 3, 4, 1, 1, 1, 5, 5, 1, 2, 3, 1, 3, 3, 1, 2, 2, 1, 2, 3, 3, 1, 1}; */
	/* int r_colors[44] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; */
	/* int c_nums[42] = {3, 1, 6, 1, 2, 3, 2, 6, 3, 7, 5, 4, 2, 3, 1, 2, 1, 1, 2, 2, 3, 5, 7, 4, 2, 1, 4, 3, 1, 1, 1, 1, 5, 3, 2, 2, 6, 1, 1, 2, 1, 1}; */
	/* int c_colors[42] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; */

	puzzle p = puzzle(w, h, &*r_sizes.begin(), &*c_sizes.begin(),
		&*r_nums.begin(), &*r_colors.begin(), &*c_nums.begin(), &*c_colors.begin());

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
	solve(&p, s);

	bool correct = check_solution(&p, s);
	if (correct)
		printf("CORRECTNESS PASSED\n");
	else
		printf("CORRECTNESS FAILED\n");

	return 0;
}
