#include <stdio.h>
#include <stdlib.h>
#include <fstream>	// File reading
#include <vector>
#include <string.h>
#include <sstream>
#include <iostream>

using namespace std;

#include "puzzle.h"
#include "solver.h"
#include "../util/CycleTimer.h"

int main(int argc, char *argv[]) {
    if(argc != 2) { printf("Expected one file name!\n"); return 0; }

    ifstream myfile;
    string line;
    bool row = true;


    myfile.open(argv[1]);
    int w = 0, h = 0;
    vector<int> r_sizes, r_nums, r_colors, c_sizes, c_nums, c_colors;
    // Export from here: http://webpbn.com/export.cgi 
    // Option: Format for Andrew Makhorin's pbnsol nonogram solver. Not for multi-color puzzles.
    while (getline(myfile,line))
    {
        if (line[0] == '*') { continue; }
        if (line[0] == '&') { row = false; continue; }
        if (row) {
            h++;
            string buf;	// Have a buffer string
            stringstream ss(line);	// Insert the string into a stream
            vector<string> tokens;	// Create vector to hold our words
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
    printf("Height: %d\n", h);
    printf("Width: %d\n",w);
    myfile.close();

	puzzle p = puzzle(w, h, &*r_sizes.begin(), &*c_sizes.begin(),
		&*r_nums.begin(), &*r_colors.begin(), &*c_nums.begin(), &*c_colors.begin());

	// Print column constraints
	printf("Column constraints:\n");
	for (int i = 0; i < p.width; i++) {
		for (int j = 0; j < p.col_sizes[i]; j++) {
			printf("%d ", p.col_constraints[i][j].num);
		}
		printf("\n");
	}
	printf("\n");
	printf("Row constraints:\n");
	// Print row constraints
	for (int i = 0; i < p.height; i++) {
		for (int j = 0; j < p.row_sizes[i]; j++) {
			printf("%d ", p.row_constraints[i][j].num);
		}
		printf("\n");
	}
	printf("\n");

	Solution s = initialize_solution(w, h);
	double start = CycleTimer::currentSeconds();
	solve(&p, s);
	double ref_time = CycleTimer::currentSeconds() - start;

	bool correct = check_solution(&p, s);

    printf("Solution:\n");
    s->print_solution();

	if (correct)
		printf("Correctness passed!\n");
	else
		printf("CORRECTNESS FAILED\n");
	printf("TIME: %.4f\n", ref_time);
    
    ofstream output; 
    output.open("output.pbm");

    int scale = 4;
    int sw = w*scale;
    int sh = h*scale;
    int8_t *d = s->data;
    int scaled[scale * scale * w* h]; 
    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            int num = d[i*w+j];
            int row = i*scale;
            int col = j*scale;
            for(int k = 0; k < scale; k++){
                for(int l = 0; l < scale; l++){
                    scaled[(row+k)*sw + col + l] = num;
                }
            }
        }
    }
    
    output << "P1\n";
    output << sw;
    output << " ";
    output << sh; 
    output << "\n"; 

    for(int i = 0; i < sh; i++){
        for(int j = 0; j < sw; j++){
            output << scaled[i*sw+j];
            output << " ";
        }
        output << "\n";
    }
    output.close();
	return 0;
}
