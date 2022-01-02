/**
 * @file TravellingSalesman.cpp
 * @author Ghazal Khalili (khalilig@mcmaster.ca)
 * @brief Solving Traveling Salesman Problem using Genetic Algorithm and Simulated Annealing
 * @version 1.1
 * @date 2022-01-01
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <cmath>
#include <numeric>
#include <iomanip>
#include "timer.hpp"
#include "matrix.hpp"
#include <random>
#include <algorithm>
#include "MyProgressbar.hpp"
#include "MyConfiguration.hpp"
#include "ChromosomeClass.hpp"
#include "ReadCSV.hpp"

using namespace std;

/**
 * @brief Overloaded binary operator << to easily print out a vector to a stream.
 * @tparam T
 * @param out The out stream.
 * @param vec The vector to print.
 * @return ostream&
 */
template <typename T>
ostream &operator<<(ostream &out, const vector<T> &vec)
{
    for (const T &i : vec)
    {
        out << i << "  ";
    }
    out << "\n";
    return out;
}

/**
 * @brief Calculates Euclidean distance.
 * @tparam T
 * @param m A matrix of points with their coordinates.
 * @return matrix<T>: Containing the length of every possible line between the given points.
 */
template <typename T>
matrix<T> Euclidean_distance(const matrix<T> &m)
{
    uint64_t row_size = m.get_rows();

    // Total number of combinations of pair of nodes to shape a line
    uint64_t combinations = (row_size - 1) * (row_size) / 2;
    vector<T> m_elems(1); // Vector of the distance matrix
    m_elems.clear();
    m_elems.reserve(combinations * 3); // Three columns (i, j, distance)
    T dis = 0;
    for (uint64_t i = 0; i < row_size; i++)
    {
        for (uint64_t j = i + 1; j < row_size; j++)
        {
            dis = sqrt((m(i, 1) - m(j, 1)) * (m(i, 1) - m(j, 1)) + (m(i, 2) - m(j, 2)) * (m(i, 2) - m(j, 2)));
            m_elems.insert(m_elems.end(), {(T)(i + 1), (T)(j + 1), dis});
        }
    }
    matrix<T> euc_distances(combinations, 3, m_elems); // The matrix of distances
    return euc_distances;
}

/**
 * @brief Exception to be thrown if the non-negativity condition of a parameter is not satisfied.
 */
class par_negative : public invalid_argument
{
public:
    par_negative() : invalid_argument("The parameter cannot be negative!\n"){};
};

/**
 * @brief Exception to be thrown if the integer condition of a parameter is not satisfied.
 */
class par_not_integer : public invalid_argument
{
public:
    par_not_integer() : invalid_argument("The parameter must be an integer!\n"){};
};

int main()
{
    /**
     * @mainpage
     * A dynamic method to solve TSP using the combination of Genetic algorithm and Simulated annealing algorithm
     */

    timer t;                 // An object helping to calculate the run time
    progressbar<double> prg; // An object to print the progress in descrete steps

    // The path of the data set file
    string datafile;
    // Defining the parameters of our algorithm (the default values are given)
    // The initial population size
    uint64_t pop_size;
    // The algorithm iterates n times
    uint64_t n_iteration;
    // Parameters needed for annealing
    uint64_t k_max; // Number of iterations
    double Temp0;   // Initial temperature
    double alpha;   // Cool-down rate

    try
    {
        // Reading the configuration file
        config configuration("config.cfg");
        datafile = configuration.get_value_string("datafile");
        double signed_v = configuration.get_value("pop_size", 50);
        if (signed_v < 0)
        {
            cout << "pop_size: ";
            throw par_negative();
        }
        if (signed_v - floor(signed_v))
        {
            cout << "pop_size: ";
            throw par_not_integer();
        }
        pop_size = (uint64_t)signed_v;
        signed_v = configuration.get_value("n_iteration", 100);
        if (signed_v < 0)
        {
            cout << "n_iteration: ";
            throw par_negative();
        }
        if (signed_v - floor(signed_v))
        {
            cout << "n_iteration:";
            throw par_not_integer();
        }
        n_iteration = (uint64_t)signed_v;
        signed_v = configuration.get_value("k_max", 100);
        if (signed_v < 0)
        {
            cout << "k_max: ";
            throw par_negative();
        }
        if (signed_v - floor(signed_v))
        {
            cout << "k_max: ";
            throw par_not_integer();
        }
        k_max = (uint64_t)signed_v;
        Temp0 = configuration.get_value("Temp0", 10);
        if (Temp0 < 0)
        {
            cout << "Temp0: ";
            throw par_negative();
        }
        alpha = configuration.get_value("alpha", 0.95);
        if (alpha < 0)
        {
            cout << "alpha: ";
            throw par_negative();
        }
    }
    catch (const exception &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }

    // Reading the data
    uint64_t Nrows, Ncols;
    try
    {
        csv<double> dataset(datafile);
        Nrows = dataset.get_NRows();
        Ncols = dataset.get_NCols();
    }
    catch (const exception &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }
    // Saving our data set into a matrix
    matrix<double> GeneData(Nrows, (Ncols + 1));
    vector<double> nodes;
    try
    {
        csv<double> dataset(datafile);
        GeneData = dataset.read_data();
        nodes = dataset.get_row_numbers();
    }
    catch (const exception &e)
    {
        cout << e.what();
        return -1;
    }
    // Updating and printing the progress
    prg.update_progress(10);
    prg.print_progress();

    // Solving TSP using the GA
    // Calculating distances (objective values) between the nodes
    matrix<double> m_distances = Euclidean_distance(GeneData);
    // Updating and printing the progress
    prg.update_progress(25);
    prg.print_progress();

    // Generating an initial population of size "pop_size"
    
    vector<chromosome<double>> population = Pop_Generator(nodes, pop_size);
    cout << "\nPrior to any optimization, ";
    double sum_ov = 0;
    for (chromosome<double> &i : population)
    {
        sum_ov += objective_value(i, m_distances);
    }
    sum_ov /= (double)pop_size;
    cout << "the mean of route length in the initial population was:\n"
         << sum_ov << "\n";

    pop_size *= 2;
    population.reserve(pop_size);
    // Updating and printing the progress
    prg.update_progress(30);
    prg.print_progress();
    double progress = ((95 - 30) / (double)n_iteration); // The main process is about 65% of the algorithm load
    uint64_t n = n_iteration;
    // Repeating n times
    while (n--)
    {
        // Regenerating another pop_size
        // (50% with cross_over, 20% mutation, 15% simulated anealling, and 15% new random)
        uint64_t pop50 = (uint64_t)round(0.5 * (double)pop_size / 2);
        if (pop50 % 2) // We need pairs of chromosomes to create children, so it should be an even number
        {
            pop50++;
        }
        uint64_t pop20 = (uint64_t)round(0.20 * (double)pop_size / 2);
        if (pop20 % 2)
        {
            pop20++;
        }
        uint64_t pop15A = (uint64_t)round(0.15 * (double)pop_size / 2);
        uint64_t pop15r = pop_size / 2 - pop20 - pop15A - pop50;
        uint64_t p = 0;       // The index for the population
        uint64_t counter = 0; // A counter of children (new offspring chromosomes)
        try
        {
            while (counter < pop50) // Cross over
            {
                vector<chromosome<double>> children = population[p].cross_over(population[p + 1]);
                if (!children[0].is_in_vector(population) && !children[1].is_in_vector(population))
                {
                    population.push_back(children[0]);
                    population.push_back(children[1]);
                    p += 2;
                    counter += 2;
                }
            }
            counter = 0;
            while (counter < pop20) // Mutation
            {
                vector<chromosome<double>> children = population[p].mutation(population[p + 1]);
                if (!children[0].is_in_vector(population) && !children[1].is_in_vector(population))
                {
                    population.push_back(children[0]);
                    population.push_back(children[1]);
                    p += 2;
                    counter += 2;
                }
            }
            counter = 0;
            while (counter < pop15A) // Simulated Annealing
            {
                chromosome<double> child = population[p].annealing(m_distances, k_max, Temp0, alpha);
                if (!child.is_in_vector(population))
                {
                    population.push_back(child);
                    p++;
                    counter++;
                }
            }
            counter = 0;
            while (counter < pop15r) // New random generation
            {
                vector<chromosome<double>> rand_population = Pop_Generator(nodes, pop15r - counter);
                for (chromosome<double> i : rand_population)
                {
                    if (!i.is_in_vector(population))
                    {
                        population.push_back(i);
                        counter++;
                    }
                }
            }

            // Selection phase
            // Calculating the objective value (total distance/cost in TSP) for each individual
            for (chromosome<double> &i : population)
            {
                i.update_fitness(1 / objective_value(i, m_distances));
            }

            if (nodes.size() <= 50)
            {
                if (n % 5 == 0) // Updating and printing the progress every 54 iterations
                {
                    prg.update_progress(30 + progress * (double)(n_iteration - n));
                    prg.print_progress();
                }
            }
            else
            { // Updating and printing the progress every iteration
                prg.update_progress(30 + progress * (double)(n_iteration - n));
                prg.print_progress();
            }
        }
        // Chatching any exception if occurred during the chromosome creation
        catch (const exception &e)
        {
            cout << "\nError: " << e.what() << '\n';
            return -1;
        }

        // Choosing 50% of the population with roulette wheel selection mechanism
        double total_objective_values = 0;
        // Calculating the chance of each chromosome to be selected
        for (const chromosome<double> &i : population)
        {
            total_objective_values += i.get_fitness();
        }
        vector<double> roulette_wheel;
        roulette_wheel.reserve(pop_size);
        double cumulative = 0;
        for (uint64_t i = 0; i < pop_size; i++)
        {
            double chance = (double)round(population[i].get_fitness() / total_objective_values * 100);
            cumulative += chance;
            roulette_wheel.push_back(cumulative);
        }
        // Selecting half of the population
        double total_chance = roulette_wheel[pop_size - 1];
        vector<chromosome<double>> selected_population;
        selected_population.reserve(pop_size / 2);
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> urd(0, total_chance);
        uint64_t found = 0;
        for (uint64_t j = 0; j < pop_size / 2; j++)
        {
            uint64_t added = 0;
            double k = urd(mt);
            for (int64_t i = (pop_size - 2); i >= 0; i--)
            {
                if (k > roulette_wheel[i])
                {
                    found = 1;
                    if (!(population[i + 1].is_in_vector(selected_population)))
                    { // The chromosome has not been added before
                        selected_population.push_back(population[i + 1]);
                        added = 1;
                        break;
                    }
                }
            }
            if (!found)
            {
                selected_population.push_back(population[0]);
                found = 1;
                added = 1;
            }
            if (added == 0)
            {
                j--;
            }
        }
        // Updating the population
        population = selected_population;
    }

    // Now choosing the best among the last population (the optimal solution)
    uint64_t max_index = find_max(population);
    // Updating and printing the progress
    prg.update_progress(100);
    prg.print_progress();
    // Printing the estimation of the optimal solution
    cout << "\nOptimal solution is estimated to be:\n"
         << population[max_index] << "with total distance of: "
         << 1 / population[max_index].get_fitness() << " units.\n";
    t.end();
    cout << "\nElapsed time: " << t.seconds() << " seconds.\n\n";
}
