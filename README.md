# Introduction

In graph theory, Travelling Salesman Problem (TSP) is a popular problem, which tries to find the optimal path. The goal is to start from a point and finish at the same starting point, crossing all the other points only once. In other words, a salesman wants to visit each and every city and come back to the original city. In TSP, the optimal solution is usually considered to be the shortest route.

TSP is an optimization problem, which tries to minimize the total covered distance. There are a variety of optimization algorithms, such as Particle Swarm, Ant Colony, Bees Colony and Genetic Algorithms, which are inspired by the behaviour of nature. 

Genetic Algorithm tries to search the solution space as much as possible, which is an excellent feature in finding the optimal solution in larger dimensions. However, it might be stuck in the local optimum. 

Simulated Annealing is another optimization algorithm inspired by the heating and controlled cooling of metal in the annealing process. The advantage of this algorithm is that it tries to scape the local optimum and reach the global optimum.

We wish to benefit from both algorithms, and we propose a dynamic method to solve TSP using the Genetic algorithm and Simulated annealing algorithm. We code the dynamic algorithm in c++.

# Genetic Algorithm (GA)
Genetic Algorithm is an effort at mimicking biology. It produces a random population of chromosomes, and each chromosome is a possible solution to the problem. The genes of the chromosome are the variables (each node in TSP).  The chromosomes cross over (recombine) and go through mutation to create children solutions (offspring chromosomes). Among all the available chromosomes, the fittest (with the shorter paths in TSP) have a higher chance to remain. The procedure iterates until reaching a reasonable estimate of the optimal solution.

# Simulated Annealing (SA)
The idea of Simulated Annealing comes from the metallurgy process, which tries to improve the physical properties of a material. The algorithm has a starting temperature and reduces the temperature by a rate of $\alpha$. During the colling, the algorithm finds a neighbour of the current state, and it accepts it if it improves the objective function. Otherwise, it accepts it by a probability value. The energy of the system gradually decreases until reaching the desired state.

# Dataset

The dataset is retrieved from [Fun TSP Challenge](https://github.com/acu192/fun-tsp-challenge). It contains the x and y coordinates of several cities. We chose to run our code on the "tiny" and "medium" datasets. We modified them a little bit and added the column names of "x" and "y". At the end of the code, the total time elapsed is printed. 

# Methodology
Our code mainly depends on a class of objects called "Chromosome". The main objective of chromosome class is to create new chromosomes. Using GA, we defined three creation methods. Also, using SA, we defined another method for modifying the chromosomes. 

- **Random creation**: It is mainly used for the initial population. The inner genes are scrambled randomly. However, the first and last nodes are set to be 1 (the constraint of the TSP of returnning to the first city). We dedicated 15% of the new population to the random creation in each iteration. This ensures that the algorithm searches the solution space as much as possible.
  
- **Crossover**: We used the single-point crossover method. In this method, two chromosome parents split from a random point, switching their genes from that point. In each child, the new part of genes from the second parent is shuffled randomly. We dedicated 50% of the new population in each iteration to this method.

- **Mutation**: We used the swap mutation method. In this method, two random genes of the offspring chromosomes are swapped. This ensures enough diversity in the population. We dedicated 20% of the new population in each iteration to this method.

- **Annealing**: We used the inversion method to find a new neighbour of the current chromosome. In this method, the genes between two randomly selected points are inversed.  In each iteration, 15% of the population is improved using the annealing method. For the parameters, we set the initial temperature to 10 and decrease it by a rate of $\alpha^{(k-1)}-1$, where $\alpha$ is 0.95, and $k$ is the current iteration. We iterated 100 times.

Finally, for the Selection phase, we employ a roulette wheel approach. It assigns higher probabilities to better solutions (fittest chromosomes).

We repeated creating and improving the population 100 times. Finally, the best chromosome is selected.

# Results
First, we should notice that the algorithms do not guarantee the global optimum, and it is an estimation of the optimum solution. Prior to any optimization, the mean of the route length was around 22 units for the tiny dataset, and 48 units for the medium dataset.

We ran the code on a system with Intel(R) Core(TM) i3-7100U CPU @ 2.40GHz   2.40 GHz processor and 8.00 GB (7.89 GB usable) of RAM. On the tiny dataset with ten cities, the result of one run was: 

```cpp
1  10  3  9  6  2  8  7  4  5  1

with total distance of: 15.1907 units.

Elapsed time: 5.41833 seconds.
```

At another run, we got:

```cpp
1  4  10  5  2  6  9  3  8  7  1  

with total distance of: 13.6214 units.

Elapsed time: 5.76986 seconds.
```
We can observe some similarities between each run. For example, in both of them sets of {7, 8}, {3, 9} ,{9, 6}, and {2, 6} appeared in an adjacent order.

Also, for the medium dataset with a new set of 100 cities, we got:
```cpp
Found the number of rows and columns:
         Rows: 100      Columns: 2

Started reading the data: **********
Reached end of the file.


[ ■■■■■                ] progress: 25%

[ ■■■■■■               ] progress: 30%

Prior to any optimization, the mean of route length in the initial population was:
48.3653

[ ■■■■■■■              ] progress: 35%

[ ■■■■■■■■             ] progress: 40%

[ ■■■■■■■■■            ] progress: 45%

[ ■■■■■■■■■■           ] progress: 50%

[ ■■■■■■■■■■■          ] progress: 55%

[ ■■■■■■■■■■■■         ] progress: 60%

[ ■■■■■■■■■■■■■        ] progress: 65%

[ ■■■■■■■■■■■■■■       ] progress: 70%

[ ■■■■■■■■■■■■■■■      ] progress: 75%

[ ■■■■■■■■■■■■■■■■     ] progress: 80%

[ ■■■■■■■■■■■■■■■■■    ] progress: 85%

[ ■■■■■■■■■■■■■■■■■■■  ] progress: 95%

[ ■■■■■■■■■■■■■■■■■■■■ ] progress: 100%

Optimal solution is estimated to be:

1  39  43  60  45  10  65  22  77  100  28  6  14  7  19  96  47  84  12  62  36  91  86  59  93  52  97  25  4  11  56  54  66  70  74  40  29  33  79  2  15  68  57  48  67  89  13  81  8  76  88  85  37  90  51  99  95  73  30  27  94  49  55 
 26  44  17  63  3  16  9  92  38  32  80  83  71  53  50  23  87  31  5  34  58  78  98  72  82  18  20  69  21  24  64  46  75  35  41  61  42  1

with total distance of: 42.3304 units.

Elapsed time: 1470.78 seconds.
```

We can see that it performs considerably slower. The reason is that there are many possible links between nodes, and trying to find them and calculating the fitness values take much more time.

# Credits
The *matrix.hpp* and *timer.hpp* are retrieved from Shoshany, B (2021) 
multithreaded-matrix (Version 1.1) [Github repository](https://github.com/bshoshany/multithreaded-matrix), and Shoshany, B [Lecture Notes for CSE 701: Foundations of Modern Scientific Programming](https://baraksh.com/CSE701/notes.php) section 7.3.3 ("Interlude: measuring performance with chrono"), respectively.
