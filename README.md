# Solving Traveling Salesman Problem using Genetic Algorithm and Simulated Annealing

- [Introduction](#introduction)
- [Genetic Algorithm (GA)](#genetic-algorithm-ga)
- [Simulated Annealing (SA)](#simulated-annealing-sa)
- [Hybrid Mechanism](#hybrid-mechanism)
- [Methodology](#methodology)
- [Usage](#usage)
- [Datasets (Sample Inputs)](#datasets-sample-inputs)
- [Sample Results](#sample-results)
- [Version history](#version-history)
- [Acknowledgment](#acknowledgment)
- [Feedback](#feedback)
- [Author and copyright](#author-and-copyright)

## Introduction

Travelling Salesman Problem (TSP) is a well-known problem in the Combinatorial Optimization field. In this problem, a salesman wants to visit a set of cities or points and then return to its original point. So, the main goal is to determine the optimal route (order of points) such that the total travelled distance is minimized, and each point is crossed only once. In graph theory, this path is called the shortest Hamiltonian Cycle. The TSP problem seems an easy-to-solve problem at first glance; however, this problem is categorized into NP-Complete problems, which means that there is no polynomial-time algorithm to exactly solve the problem. Therefore, when the dimensionality of the problem grows, the execution time increases significantly. With this in mind, scholars and scientists use randomized optimization algorithms to obtain good solutions for TSP in a reasonable time.

There are some famous randomized algorithms called Meta-heuristics, some of which are stimulated by observations from nature. Ants colony optimization algorithm, Artificial bee colony algorithm, Swarm intelligence, Genetic algorithm, Simulated annealing, and Tabu search are some examples of meta-heuristic algorithms.

Genetic Algorithm tries to search the solution space as much as possible, which is an excellent feature in finding the optimal solution in larger dimensions. However, it might be stuck in the local optimum.

Simulated Annealing is another optimization algorithm inspired by the heating and controlled cooling of metal in the annealing process. The advantage of this algorithm is that it tries to scape the local optimum and reach the global optimum.

We wish to benefit from both algorithms, and we propose a dynamic method to solve TSP using the combination of Genetic algorithm and Simulated annealing algorithm. We code the dynamic algorithm in c++.

## Genetic Algorithm (GA)

Genetic Algorithm is an effort at mimicking biology and is inspired by chromosomes' behavior. Every chromosome comprises several genes defining the characteristics of the entire chromosome. More robust characteristics are more likely to be inherited by the next generation.

To solve TSP, we consider each chromosome as a possible solution to the problem, and the genes are the points (nodes of a graph). In the initial stage, a random population of chromosomes is created. Then, the population undergoes two main processes:

1. Crossover (recombine): Pairs of parent chromosomes create the next generation of chromosomes, known as children or offspring chromosomes. There are plenty of approaches for the crossover, among which we use the single-point method.
2. Mutation: A specific portion of offspring chromosomes go through mutation, and some of their characteristics change. Likewise, there are many methods for the mutation, among which we take advantage of the swap mutation approach.

After these processes, the fittest chromosomes remain to produce the next generation. For this purpose, we use the roulette wheel. In this method, those chromosomes with a smaller total distance have a higher chance of being selected. Afterwards, this procedure iterates until reaching a reasonable estimate of the optimal solution.

## Simulated Annealing (SA)

The idea of Simulated Annealing comes from the metallurgy process, which tries to improve the physical properties of a material. Annealing is a heat treatment process whereby the hardness of materials decreases. Therefore, the material becomes more workable. Likewise, Simulated Annealing is an artificial annealing process for mathematical optimization. Although the main idea of this algorithm seems simple, simulated annealing has shown appropriate performance for solving NP-Complete problems.

The algorithm has a starting temperature and reduces the temperature by a rate of $\alpha$. During the cooling, the algorithm finds a neighbor of the current state, and it accepts it if it improves the objective function. Otherwise, it accepts it by a probability value.

In other words, the current solution undergoes a change. Afterwards, the current solution will be replaced with the new solution if the new one is better than the current solution. Otherwise, the newly found solution is accepted probabilistically. This procedure occurs iteratively. In the initial steps, the probability of accepting a worse solution is high. However, this probability decreases as the algorithm goes on, resembling the annealing process for materials. The main goal of this approach is to escape local optimums and find better solutions.

## Hybrid Mechanism

For solving NP-Complete problems, it is customary to mix meta-heuristic algorithms to get better performance. Therefore, we take the same strategy and embed the simulated annealing inside the genetic algorithm. In other words, after the cross-over and mutation processes, we dedicated a particular portion of the chromosome population to the simulated annealing. The selected chromosomes, one by one, will be considered as the starting point for the simulated annealing algorithm. Then, the improved solutions will be added to the main population of chromosomes, and the genetic algorithm goes on.

## Methodology

Our code mainly depends on a class of objects called `Chromosome`. The main objective of chromosome class is to create new chromosomes. By default, the algorithm generates 50 random chromosomes as the initial population. Using GA, we defined three creation methods; and using SA, we defined another method for modifying the chromosomes.

- **Random creation**: It is mainly used for the initial population. The inner genes are scrambled randomly. However, the first and last nodes are set to be 1 (the constraint of the TSP of returning to the first city). We dedicated 15% of the new population to the random creation in each iteration. This ensures that the algorithm searches the solution space as much as possible.
  
- **Crossover**: We used the single-point crossover method. In this method, two parent chromosomes split from a random point, switching their genes from that point. In each child, the new part of genes from the second parent is shuffled randomly. We dedicated 50% of the new population in each iteration to this method.

- **Mutation**: We used the swap mutation method. In this method, two random genes of the offspring chromosomes are swapped. This ensures enough diversity in the population. We dedicated 20% of the new population in each iteration to this method.

- **Annealing**: We used the inversion method to find a new neighbor of the current chromosome. In this method, the genes between two randomly selected points are inverted. In each iteration, 15% of the population is improved using the annealing method. For the parameters, we set the default initial temperature to be 10 and decrease it by a rate of $\alpha^{(k-1)}-1$, where $\alpha$ is 0.95, and $k$ is the current iteration. By default, the algorithm iterates 100 times.

Finally, for the selection phase, we employ a roulette wheel approach. It assigns higher probabilities to better solutions (fittest chromosomes).

The algorithm repeats creating and improving the population 100 times, by default. Finally, the best chromosome is selected.

## Usage

Prior to run the `TravellingSalesman.cpp`, you need to include the following header files:

- `ChromosomeClass.hpp`: Contains the `chromosome` class, with member functions to generate and improve the chromosome population.
- `MyConfiguration.hpp`: Contains the `config` class for reading the configuration file. The parameters are defined in the format of `key = value; comments`. If a parameter is not defined, the default value is used.
- `ReadCSV.hpp`: Contains the `csv` class for reading the input dataset file in csv format.
- `MyProgressbar.hpp`: Contains the `progressbar` class for printing the progress bar in discrete steps.
- `matrix.hpp`: Contains the `matrix` class, with member functions for saving and accessing the data in a matrix format.
- `timer.hpp`: Contains the `timer` class for calculating the elapsed time in running the algorithm.

## Datasets (Sample Inputs)

The dataset is retrieved from [Fun TSP Challenge](https://github.com/acu192/fun-tsp-challenge). It contains the x and y coordinates of several cities. We chose to run our code on the 'tiny' and 'medium' datasets. We modified them a little bit and added the column names of 'x' and 'y'. At the end of the code, the total time elapsed is printed.

## Sample Results

We should notice that the algorithms do not guarantee the global optimum, and it is an estimation of the optimum solution. Prior to any optimization, the mean of the route length was around 22 units for the tiny dataset, and 48 units for the medium dataset. The input dataset can be altered by changing the configuration file.

We ran the code on a system with Intel(R) Core(TM) i3-7100U CPU @ 2.40GHz 2.40 GHz processor and 8.00 GB (7.89 GB usable) of RAM, using GCC compiler v11.2.0 on Windows 10 build 19044.1415. To satisfy the C++20 requirement, the compiler must have the `-std=c++20` flag. On the tiny dataset with ten cities, the result of one run was:

```none
Configuration file is received

Data file is successfully received

Found the number of rows and columns:
         Rows: 10       Columns: 2

Started reading the data: **********
Reached end of the file.
 + All the rows are received successfully.


[ *****                                              ] progress: 10%

[ ************                                       ] progress: 25%

Prior to any optimization, the mean of route length in the initial population was:
22.8937

[ ***************                                    ] progress: 30%

[ ****************                                   ] progress: 33.25%

[ ******************                                 ] progress: 36.5%

[ *******************                                ] progress: 39.75%

[ *********************                              ] progress: 43%

[ ***********************                            ] progress: 46.25%

[ ************************                           ] progress: 49.5%

[ **************************                         ] progress: 52.75%

[ ****************************                       ] progress: 56%

[ *****************************                      ] progress: 59.25%

[ *******************************                    ] progress: 62.5%

[ ********************************                   ] progress: 65.75%

[ **********************************                 ] progress: 69%

[ ************************************               ] progress: 72.25%

[ *************************************              ] progress: 75.5%

[ ***************************************            ] progress: 78.75%

[ *****************************************          ] progress: 82%

[ ******************************************         ] progress: 85.25%

[ ********************************************       ] progress: 88.5%

[ *********************************************      ] progress: 91.75%

[ ***********************************************    ] progress: 95%

[ ************************************************** ] progress: 100%

Optimal solution is estimated to be:

1  4  5  2  6  9  3  8  7  10  1

with total distance of: 12.7852 units.

Elapsed time: 6.90059 seconds.
```

At another run, we got:

```none
Optimal solution is estimated to be:

1  5  2  6  9  8  7  3  10  4  1  

with total distance of: 14.4916 units.

Elapsed time: 7.12784 seconds.
```

We can observe some similarities between each run. For example, in both of them sets of {1, 4}, {2, 5}, {2, 6}, {6, 9}, and {7, 8} appeared in an adjacent order.

Also, for the medium dataset with a new set of 100 cities, we got:

```none

Prior to any optimization, the mean of route length in the initial population was:
47.9733

[ ************************************************** ] progress: 100%

Optimal solution is estimated to be:

1  28  33  39  68  6  20  13  30  9  81  57  2  31  36  35  76  92  27  63  23  80  34  37  49  5  99  17  38  11  26  95  
61  97  21  66  29  64  87  71  7  46  79  96  53  93  82  75  67  47  58  88  72  54  12  59  42  41  19  62  69  70  78  
40  73  52  45  91  56  44  100  15  83  94  43  50  22  85  77  74  14  84  98  65  8  10  16  60  51  86  18  3  55  48  
90  32  89  24  4  25  1

with total distance of: 41.3113 units.

Elapsed time: 1418.13 seconds.
```

We can see that it performs considerably slower. The reason is that there are many possible links between nodes, and trying to find them and calculating the fitness values take much more time.

## Version history

- Version 1.1 (2022-01-01)
  - Some cosmetic changes.
  - Used updated versions of `MyConfiguration.hpp` and `ReadCSV.hpp`, and edited reading the files.
  - The `README.md` documentation file is updated.
  - [Repository](https://github.com/Ghzlikli/TSP-Hybrid_Optimization) renamed.
- Version 1.0 (2021-12-23)
  - Initial release.

## Acknowledgment

The `matrix.hpp` and `timer.hpp` are retrieved from Shoshany, Barak [Lecture Notes for CSE 701: Foundations of Modern Scientific Programming](https://baraksh.com/CSE701/notes.php) section 7.1.3 ('Class template example: the matrix class again'), and section 7.3.3 ('Interlude: measuring performance with chrono'), respectively.

The `MyConfiguration.hpp` and `ReadCSV.hpp` are retrieved from [Configuration](https://github.com/Ghzlikli/Configuration) and [Reading-CSV](https://github.com/Ghzlikli/Reading-CSV) repositories, respectively.

## Feedback

Any feedback or request is highly appreciated through [open a new issue](https://github.com/Ghzlikli/TSP-Hybrid_Optimization/issues) or [email](mailto:khalili.ghazal.97@gmail.com).

## Author and copyright

Copyright (c) 2021 [Ghazal Khalili](mailto:khalili.ghazal.97@gmail.com).

If you use this class template in your code, please acknowledge the author and provide a link to the [GitHub repository](https://github.com/Ghzlikli/TSP-Hybrid_Optimization). 

**Thank you!**
