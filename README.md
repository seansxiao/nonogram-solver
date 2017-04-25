## Parallel Nonogram Solver

Team members: Sean Xiao (ssx), Nicholas Mu (nhmu)

### Summary:
We will implement a nonogram solver on the GHC machines and explore the benefit of parallelizing the solver algorithm.

### Background:
[Nonograms](https://en.wikipedia.org/wiki/Nonogram) are picture logic puzzles, usually on a 2-D grid, that when solved, produces a black/white pixel image, although there are also versions that produce images with a wider color palette.

Ex:

![nonogram_example](https://upload.wikimedia.org/wikipedia/commons/d/d0/Nonogram2.jpg)

The numbers to the side of a particular row or above a particular column specify the number of contiguous blocks of color there are in that row or column, with contiguous blocks in the same row separated by at least one blank square.  In this example, the fourth row reads “2 2”, specifying that in the fourth row there are two contiguous blocks of black pixels with length 2 each, separated by (in this case) two white pixels.  Usually nonogram puzzles have a unique solution, and have no theoretical bounds on size.

There are other nonogram solver implementations online, but they are mostly sequential implementations.  We think that developing a parallel implementation (for example, working on several different sections of the image at the same time) may improve the runtime of finding a solution, particularly on larger image sizes.

### The Challenge:
There are many different algorithms and approaches to nonogram solving, and at this point our team is not too familiar with how these different algorithms work.  However, this means that there are many different ways to think about the problem apart from a simple naive or brute force solution.  It will be worthwhile for us to research these different approaches and analyze and identify points of parallelization for each in order for us to arrive at a solver that makes good use of parallel computing resources.

### Resources:
There are multiple online solvers that we plan on using to test our correctness against such as http://a.teall.info/nonogram/. We also plan on implementing some ideas discussed in this paper: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.177.76&rep=rep1&type=pdf.

### Goals and Deliverables:
Our goal is to identify or come up with a highly parallelizable nonogram solving algorithm and then to implement it in C++ on the 16-threaded GHC machines.  We plan to achieve notable speedup over the sequential implementation of the same algorithm, especially on larger input sizes, although linear speedup seems unlikely due to some of the dependencies involved in solving the puzzles as well as algorithmic specifics.  If possible, we would like our implementation to not just achieve good speedup relative to a set baseline sequential algorithm, but also be fast relative to all of the other algorithms that exist.

The project can be demo-ed in person by feeding it puzzles and then running the parallel and sequential versions of the solver and comparing the runtimes.  We will produce speedup graphs for various thread counts to demonstrate how much speedup is achieved on various puzzle inputs.

### Platform Choice:
We plan on using the GHC machines. These provide 16 execution contexts and 8-wide vector instructions which provide enough resources to test and analyze our performance on.

### Preliminary Schedule:
April 9-15: Manually practice and gain familiarity with solving nonograms. Implement correct sequential algorithms.

April 16-22: Identify areas of parallelism and finish implementing naive parallel solution.

April 23-29: Optimize first parallel solution and begin implementing different parallel solutions.

April 30 - May 6: Finish other parallel solutions.

May 7 - May 12: Analyze results and finish final write-up.

## Checkpoint (April 25)
Writing a nonogram solver turned out to be harder than we expected. So far we have almost completed the logical solving portion of the program, and it solves a large number of puzzles that only require line-solving techniques (i.e. no guessing). However, in order to solve the rest of the puzzles we have to implement a sort of DFS, which we have yet to do. So far the code is entirely sequential, and we plan on finishing the line solver and DFS before attempting to parallelize the algorithm.

With regard to the goals and deliverables, we have yet to start working on parallelization, but we still expect to be able to accomplish decent speedup. However, because writing the solver was harder than anticipated, it will probably not be as fast as some of the fastest nonogram-solving algorithms that exist.

We will have speedup graphs to show for the parallelism competition, but we will also be able to solve different puzzles on the spot and produce their resulting images. We will be able to run both the sequential and parallel versions side by side in order to compare.

### Revised Schedule:
April 25: Finish most of the line-solver

April 26 - 28: Complete line-solving algorithm

April 29 - 30: Finish implementing DFS portion of algorithm

May 1 - May 2: Optimize sequential code

May 3 - May 5: Identify areas of parallelization and begin implementing parallel code

May 6 - May 8: Finish parallelization implementation and start optimizing

May 9 - May 12: Analyze results and finish final write-up.
