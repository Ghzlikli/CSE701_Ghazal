/**
 * @file Projet 3.cpp
 * @author Ghazal Khalili (khalilig@mcmaster.ca)
 * @brief Solving Traveling Salesman Problem using Genetic Algorithm
 * @version 0.1
 * @date 2021-12-16
 *
 * @copyright Copyright (c) 2021
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

using namespace std;

/**
 * @brief Class of Progressbar
 * to print the progress in discrete steps
 * @tparam T type which is usually unsigned integer
 */
template <typename T>
class progressbar
{
public:
    /**
     * @brief Constructs a new progressbar object.
     * with progress equals to 0
     */
    progressbar();

    /**
     * @brief Constructs a new progressbar object with a given progress.
     * @param _progress initial progress
     */
    progressbar(T const &);

    /**
     * @brief Exception to be thrown if the given progress is not between 0 and 100.
     */
    class invalid_progress : public invalid_argument
    {
    public:
        invalid_progress() : invalid_argument("Progress percentage should be between 0 and 100!"){};
    };

    /**
     * @brief Updates the current progress.
     * @param _progress the new progress
     */
    void update_progress(T const &);

    /**
     * @brief Prints the progress bar.
     */
    void print_progress();

private:
    // the current progress
    T progress = 0;
};

// ==============
// Implementation
// ==============

template <typename T>
inline progressbar<T>::progressbar() : progress(0) {}

template <typename T>
inline progressbar<T>::progressbar(T const &_progress) : progress(_progress)
{
    if (progress < 0 || progress > 100)
    {
        throw typename progressbar<T>::invalid_progress();
    }
}

template <typename T>
inline void progressbar<T>::update_progress(T const &_progress)
{
    if (_progress < 0 || _progress > 100)
    {
        throw typename progressbar<T>::invalid_progress();
    }
    progress = _progress;
}

template <typename T>
inline void progressbar<T>::print_progress()
{
    // each 5% has a square or space in the progress bar
    uint64_t filled_sq = (uint64_t)(progress / 5);
    uint64_t spaces = 20 - filled_sq;
    cout << "\n[ ";
    while (filled_sq--)
    {
        cout << "■";
    }
    while (spaces--)
    {
        cout << " ";
    }
    cout << " ] progress: " << progress << "%\n";
}

// =========================
// End of Progressbar Class
// =========================

/**
 * @brief Overloaded binary operator << to easily print out a vector to a stream.
 * @tparam T
 * @param out the out stream
 * @param vec the vector to print
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
 * @brief Reads a row of the dataset in csv format
 * (columns are separated by comma).
 * @param s the row of the data set
 * @param NCols number of columns
 * @return vector<double> the read row
 */
vector<double> read_D_rows(const string &s, const uint64_t &NCols)
{
    vector<double> v; // the vector for saving the data
    uint64_t j = 0;   // number of columns read for each row
    istringstream string_stream(s);
    string ss;
    try
    {
        // extracting values of the columns from each row
        while (getline(string_stream, ss, ','))
        {
            // checking whether it is a number
            vector<char> cstr(ss.c_str(), ss.c_str() + ss.size() + 1); // breaking down the string to characters
            uint64_t dot = 0;
            uint64_t size = cstr.size() - 1;
            for (uint64_t i = 0; i < size; i++)
            {
                if (cstr[i] == '.') // the number is allowed to have a decimal point
                {
                    if (i == size - 1) // we should have at least one digit after the decimal point
                        throw invalid_argument("Expected a number!\n");
                    dot++;
                    continue;
                }
                if (i == 0 && cstr[i] == '-') // in case of a negative number, skips the minus sign
                {
                    continue;
                }
                // the other characters should be only digits, and no more than one decimal point is allowed
                if (isdigit(cstr[i]) == 0 || dot > 1)
                {
                    throw invalid_argument("Expected a number!\n");
                }
            }
            // adding the number to our vector
            v.push_back(stod(ss));
            j++; // updating the size
        }
    }
    catch (const invalid_argument &e)
    {
        // checking whether there is not a number
        throw invalid_argument("Expected a number!\n");
    }
    catch (const out_of_range &e)
    {
        throw out_of_range("Number is out of range!\n");
    }

    // checking whether the the number of read values is correct (equals the number of columns)
    if (j < NCols)
    {
        throw invalid_argument("Expected more columns!\n");
    }
    else if (j > NCols)
    {
        throw invalid_argument("Expected less columns!\n");
    }
    return v;
}

/**
 * @brief Calculates Euclidean distance.
 * @tparam T
 * @param m a matrix of points with coordinates
 * @return matrix<T> containing the length of every possible line between the given points
 */
template <typename T>
matrix<T> Euclidean_distance(const matrix<T> &m)
{
    uint64_t row_size = m.get_rows();

    // total number of combinations of pair of nodes to shape a line
    uint64_t combinations = (row_size - 1) * (row_size) / 2;
    vector<T> m_elems(1); // vector of the distance matrix
    m_elems.clear();
    m_elems.reserve(combinations * 3); // three columns (i, j, distance)
    T dis = 0;
    for (uint64_t i = 0; i < row_size; i++)
    {
        for (uint64_t j = i + 1; j < row_size; j++)
        {
            dis = sqrt((m(i, 1) - m(j, 1)) * (m(i, 1) - m(j, 1)) + (m(i, 2) - m(j, 2)) * (m(i, 2) - m(j, 2)));
            m_elems.insert(m_elems.end(), {(T)(i + 1), (T)(j + 1), dis});
        }
    }
    matrix<T> euc_distances(combinations, 3, m_elems); // the matrix of distances
    return euc_distances;
}

/**
 * @brief Class of Chromosomes.
 * @param genes chromosomes are created by genes.
 * @param fitness proportional to objective value; inspired by the fact that
 * in nature, the fittest individuals have a higher chance of getting food and mating.
 * @tparam T
 */
template <typename T>
class chromosome
{

public:
    /**
     * @brief Constructor to create a null chromosome.
     */
    chromosome();

    /**
     * @brief Constructor to create a null chromosome of the given size.
     * @param _size the size of the new chromosome
     */
    chromosome(const uint64_t &);

    /**
     * @brief Constructor to create a new chromosome randomly.
     * @param _genes a vector of genes; the size is inferred automatically
     */
    chromosome(const vector<T> &);

    /**
     * @brief Returns size of the chromosome.
     * @return uint64_t size of the chromosome
     */
    uint64_t size() const;

    /**
     * @brief Overloaded operator [] to access genes elements.
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n the element number
     * @return T&
     */
    T &operator[](const uint64_t &);

    /**
     * @brief Overloaded operator [] to access genes elements.
     * The indices start from 0, similar to a vector;
     * Second version: does not allow modification of the element.
     * @param n the element number
     * @return T&
     */
    const T &operator[](const uint64_t &) const;

    /**
     * @brief Access genes elements with range checking (throws out_of_range via vector::at).
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n the element number
     * @return T&
     */
    T &at(const uint64_t &);

    /**
     * @brief Access genes elements with range checking (throws out_of_range via vector::at).
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n the element number
     * @return T&
     */
    const T &at(const uint64_t &) const;

    /**
     * @brief Returns a read/write iterator that points to the first element in the genes vector.
     * Iteration is done in ordinary element order.
     * @return vector<T>::iterator
     */
    vector<T>::iterator begin();

    /**
     * @brief Returns a read/write iterator that points to one past the last element in the genes vector.
     * Iteration is done in ordinary element order.
     * @return vector<T>::iterator
     */
    vector<T>::iterator end();

    /**
     * @brief Gets the genes.
     * @return vector<T> of genes
     */
    vector<T> get_genes() const;

    /**
     * @brief Updates the fitness value of the chromosome.
     * @param fit_value the given/new fitness value
     */
    void update_fitness(const double &);

    /**
     * @brief Gets the fitness value of the chromosome.
     * @return double
     */
    double get_fitness() const;

    /**
     * @brief Copies a chromosome.
     * @return chromosome<T>
     */
    chromosome<T> ch_copy();

    /**
     * @brief Genereate a specific number of random values.
     * It is needed for creating a new chromosome with genes in a random order
     * @param n number of random values needed
     * @return vector<double> the vector of n radom values
     */
    vector<double> assign_rand(const uint64_t &);

    /**
     * @brief Exception to be thrown if two chromosomes of different sizes are cross overed.
     */
    class incompatible_sizes : public invalid_argument
    {
    public:
        incompatible_sizes() : invalid_argument("Cannot continue due to different dimensions!"){};
    };

    /**
     * @brief Exception to be thrown if the distance matrix does not
     *  contain the line between the genes of the chromosome.
     */
    class no_link_found : public invalid_argument
    {
    public:
        no_link_found() : invalid_argument("The distance matrix does not contain the required link!"){};
    };

    /**
     * @brief The single point crossover.
     * @param ch the second parent
     * @return vector<chromosome<T>> two new children (offspring chromosomes)
     */
    vector<chromosome<T>> cross_over(chromosome<T> &);

    /**
     * @brief The Swap Mutation:
     * swaps two genes of the children after being thorough a cross-over.
     * @param ch the second parent
     * @return vector<chromosome<T>> two mutated children (offspring chromosomes)
     */
    vector<chromosome<T>> mutation(chromosome<T> &);

    /**
     * @brief Checks whether the chromosome is in a vector of chromosomes.
     * @param v_ch the vector of chromosomes
     * @return TRUE if it is in the container; FALSE otherwise
     */
    bool is_in_vector(const vector<chromosome<T>> &);

    /**
     * @brief Finds a chromosome from the neighborhood.
     * It can be used for Simulated Anealling
     * @return * chromosome<T> the new chromosome
     */
    chromosome<T> neighborhood();

    /**
     * @brief Runs a Simulated Annealing optimization on a given chromosome.
     * @param m_distances the matrix of distances, used for calculating the objective value
     * @return chromosome<T> the final found chromosome
     */
    chromosome<T> annealing(matrix<double> const &);

private:
    // size of the chromosome
    uint64_t ch_size = 1;
    // a vector containing the genes (nodes)
    vector<T> genes;
    /* the fitness of a chromosome to contribute to the next generation
    it is proportional to the objective value */
    double fitness = 0;
};

/**
 * @brief Overloaded binary operator << to easily print out a chromosome object to a stream.
 * @tparam T
 * @param out the out stream
 * @param ch the chromosome to print
 * @return ostream&
 */
template <typename T>
ostream &operator<<(ostream &, const chromosome<T> &);

/**
 * @brief Generates a population of new random chromosomes.
 * @tparam T
 * @param genes the genes creating the chromosomes
 * @param T_number the total number of individuals
 * @return vector<chromosome<T>>, the population (a vector of chromosomes)
 */
template <typename T>
vector<chromosome<T>> Pop_Generator(vector<T> const &, const uint64_t &);

/**
 * @brief Calculates the total cost of each solution (individual);
 * In TSP, the total cost is the total distances covered.
 * @tparam T
 * @param ch the solution (chromosome)
 * @param cost the matrix containing the cost of different lines
 * @return T
 */
template <typename T>
T objective_value(chromosome<T> &, matrix<T> const &);

/**
 * @brief Finds the individual with the maximum value of fitness
 * @tparam T
 * @param pop the population to search in
 * @return uint64_t the index of the individual with the maximum value
 */
template <typename T>
uint64_t find_max(vector<chromosome<T>> const &);

// ==============
// Implementation
// ==============

template <typename T>
inline chromosome<T>::chromosome() : genes(1)
{
}

template <typename T>
inline chromosome<T>::chromosome(const uint64_t &_size) : ch_size(_size + 1), genes(_size + 1) {}

template <typename T>
chromosome<T>::chromosome(const vector<T> &_genes) : ch_size(_genes.size() + 1)
{
    genes = vector<T>(ch_size);
    // the start and end points have to be the fixed
    genes[0] = 1;
    genes[ch_size - 1] = 1;
    vector<T> internal_nodes(ch_size - 2); // other points
    copy(_genes.begin() + 1, _genes.end(), internal_nodes.begin());
    vector<double> rand_g(ch_size - 2);
    // assign a random value for each gene
    rand_g = chromosome::assign_rand(ch_size - 2);

    uint64_t max_indx = 0; // index of the maximum value
    double swap = 0;
    T swap2 = 0;

    // in descending oreder
    for (uint64_t i = 0; i < ch_size - 3; i++)
    {
        max_indx = i;
        for (uint64_t j = i; j < ch_size - 2; j++)
        {
            if (rand_g[j] > rand_g[max_indx])
                max_indx = j;
        }
        swap = rand_g[i];
        rand_g[i] = rand_g[max_indx];
        rand_g[max_indx] = swap;

        swap2 = internal_nodes[i];
        internal_nodes[i] = internal_nodes[max_indx];
        internal_nodes[max_indx] = swap2;
    }
    copy(internal_nodes.begin(), internal_nodes.end(), genes.begin() + 1);
}

template <typename T>
inline uint64_t chromosome<T>::size() const
{
    return ch_size;
}

template <typename T>
inline T &chromosome<T>::operator[](const uint64_t &n)
{
    return genes[n];
}

template <typename T>
const inline T &chromosome<T>::operator[](const uint64_t &n) const
{
    return genes[n];
}

template <typename T>
inline T &chromosome<T>::at(const uint64_t &n)
{
    return genes.at(n);
}

template <typename T>
const inline T &chromosome<T>::at(const uint64_t &n) const
{
    return genes.at(n);
}

template <typename T>
inline vector<T> chromosome<T>::get_genes() const
{
    return genes;
}

template <typename T>
inline vector<T>::iterator chromosome<T>::begin()
{
    return genes.begin();
}

template <typename T>
inline vector<T>::iterator chromosome<T>::end()
{
    return genes.end();
}

template <typename T>
inline void chromosome<T>::update_fitness(const double &fit_value)
{
    fitness = fit_value;
}

template <typename T>
inline double chromosome<T>::get_fitness() const
{
    return fitness;
}

template <typename T>
inline chromosome<T> chromosome<T>::ch_copy()
{
    chromosome<T> new_chr(ch_size - 1);
    copy(genes.begin(), genes.end(), new_chr.begin());
    return new_chr;
}

template <typename T>
vector<double> chromosome<T>::assign_rand(const uint64_t &n)
{
    vector<double> rand_gv(n);
    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> nd(0, 1);
    for (double &i : rand_gv)
        i = nd(mt);

    return rand_gv;
}

template <typename T>
vector<chromosome<T>> chromosome<T>::cross_over(chromosome<T> &ch)
{

    uint64_t parent_size = ch.size();
    // check for equal sizes
    if (ch_size != parent_size)
    {
        throw typename chromosome<T>::incompatible_sizes();
    }

    // randomly select a point to cross over
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2);
    uint64_t k = uid(mt);

    // child 1
    chromosome<T> ch1 = ch.ch_copy();
    for (uint64_t i = k; i < ch_size - 1; i++)
    {
        for (uint64_t j = 1; j < ch_size - 1; j++)
        {
            if (ch1[j] == genes[i])
            {
                remove(ch1.begin(), ch1.end(), ch1[j]);
                break;
            }
        }
    }
    // the second part of the first parent is going to be randomly placed into the child
    vector<T> deleted_part(ch_size - k);
    copy(genes.begin() + k, genes.end(), deleted_part.begin());
    vector<double> rand_g;
    // assign a random value
    rand_g = assign_rand(deleted_part.size() - 1); // except the last 1 (ending point is fixed)
    uint64_t max_indx = 0;                         // index of the maximum value
    double swap = 0;
    T swap2 = 0;

    // in descending oreder
    for (uint64_t i = 0; i < deleted_part.size() - 1; i++)
    {
        max_indx = i;
        for (uint64_t j = i; j < deleted_part.size() - 1; j++)
        {
            if (rand_g[j] > rand_g[max_indx])
                max_indx = j;
        }
        swap = rand_g[i];
        rand_g[i] = rand_g[max_indx];
        rand_g[max_indx] = swap;

        swap2 = deleted_part[i];
        deleted_part[i] = deleted_part[max_indx];
        deleted_part[max_indx] = swap2;
    }
    copy(deleted_part.begin(), deleted_part.end(), ch1.begin() + k);

    vector<chromosome<T>> children;
    children.push_back(ch1);

    // child 2
    // copying the parents again
    chromosome<T> ch2(ch_size - 1);
    copy(genes.begin(), genes.end(), ch2.begin());
    for (uint64_t i = k; i < ch_size - 1; i++)
    {
        for (uint64_t j = 1; j < ch_size - 1; j++)
        {

            if (ch2[j] == ch[i])
            {
                remove(ch2.begin(), ch2.end(), ch2[j]);
                break;
            }
        }
    }
    // the second part of the second parent is going to be randomly placed into the child
    copy(ch.begin() + k, ch.end(), deleted_part.begin());
    // assign a random value
    rand_g = assign_rand(deleted_part.size() - 1); // except the last 1 (ending point is fixed)

    // in descending oreder
    for (uint64_t i = 0; i < deleted_part.size() - 1; i++)
    {
        max_indx = i;
        for (uint64_t j = i; j < deleted_part.size() - 1; j++)
        {
            if (rand_g[j] > rand_g[max_indx])
                max_indx = j;
        }
        swap = rand_g[i];
        rand_g[i] = rand_g[max_indx];
        rand_g[max_indx] = swap;

        swap2 = deleted_part[i];
        deleted_part[i] = deleted_part[max_indx];
        deleted_part[max_indx] = swap2;
    }
    copy(deleted_part.begin(), deleted_part.end(), ch2.begin() + k);
    children.push_back(ch2);
    return children;
}

template <typename T>
vector<chromosome<T>> chromosome<T>::mutation(chromosome<T> &ch)
{
    // copying the first parent again
    chromosome<T> parent(ch_size - 1);
    copy(genes.begin(), genes.end(), parent.begin());
    vector<chromosome<T>> children = parent.cross_over(ch);
    // randomly select two points to swap
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2); // swap except the first and last one
    uint64_t k1 = uid(mt);
    uint64_t k2 = uid(mt);

    for (uint64_t i = 0; i < 2; i++)
    {
        swap(children[i][k1], children[i][k2]);
    }
    return children;
}

template <typename T>
bool chromosome<T>::is_in_vector(const vector<chromosome<T>> &v_ch)
{
    for (const chromosome<T> &i : v_ch)
    {
        if (i.genes == genes)
        {
            return true;
        }
    }
    return false;
}

template <typename T>
chromosome<T> chromosome<T>::annealing(matrix<double> const &m_distances)
{
    uint64_t k_max = 100; // number of iterations
    double Temp0 = 10;    // initial temperature
    double Temp;
    double alpha = 0.95;
    // copying the current chromosome
    chromosome<T> ch0(ch_size - 1);
    copy(genes.begin(), genes.end(), ch0.begin());
    chromosome<T> new_ch(ch_size - 1);
    for (uint64_t k = 0; k < k_max; k++)
    {
        new_ch = ch0.neighborhood();
        double f_new = objective_value(new_ch, m_distances);
        double f0 = objective_value(ch0, m_distances);
        if (f_new < f0)
        {
            ch0 = new_ch.ch_copy(); // we accept it
        }
        else
        {
            double P = exp(-(f_new - f0) / Temp);
            // deciding on whether to accept it or not
            random_device rd;
            mt19937 mt(rd());
            uniform_real_distribution<double> urd(0, 1);
            double r = urd(mt);
            if (P >= r)
            {
                ch0 = new_ch.ch_copy(); // we accept it
            }
        }
        Temp = Temp0 * pow(alpha, k); // cooling the temperature
    }
    return ch0;
}

template <typename T>
chromosome<T> chromosome<T>::neighborhood()
{
    // copying the current chromosome
    chromosome<T> new_ch(ch_size - 1);
    copy(genes.begin(), genes.end(), new_ch.begin());
    // randomly select two points to swap
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2); // change except the first and last points
    uint64_t k1 = uid(mt);
    uint64_t k2 = uid(mt);
    reverse(new_ch.begin() + k1, new_ch.begin() + k2);
    return new_ch;
}

template <typename T>
ostream &operator<<(ostream &out, const chromosome<T> &ch)
{
    out << "\n"
        << ch.get_genes() << "\n";
    return out;
}

template <typename T>
vector<chromosome<T>> Pop_Generator(vector<T> const &genes, const uint64_t &T_number)
{
    vector<chromosome<T>> ch_Generation(T_number, chromosome<T>());
    for (uint64_t i = 0; i < T_number; i++)
    {
        chromosome<T> new_ch(genes);
        // checking whether it is already in the population or not
        if (!new_ch.is_in_vector(ch_Generation))
            ch_Generation[i] = new_ch;
        else
            i--; // repeat to get a new one
    }
    return ch_Generation;
}

template <typename T>
T objective_value(chromosome<T> &ch, matrix<T> const &cost)
{
    uint64_t size_ch = ch.size();
    uint64_t size_m = cost.get_rows();
    T fitness = 0;
    uint64_t count = 0; // counts number of lines added
    for (uint64_t i = 0; i < (size_ch - 1); i++)
    {
        T gene1 = ch[i];
        T gene2 = ch[i + 1];
        for (uint64_t j = 0; j < size_m; j++)
        {
            T node1 = cost(j, 0);
            T node2 = cost(j, 1);
            if ((gene1 == node1 && gene2 == node2) || (gene1 == node2 && gene2 == node1))
            {
                fitness += cost(j, 2);
                count++; // found another link
            }
        }
    }
    // checks if found all the required lines
    if (count != size_ch - 1)
    {
        throw typename chromosome<T>::no_link_found();
    }
    return fitness;
}

template <typename T>
uint64_t find_max(vector<chromosome<T>> const &pop)
{
    T max = pop[0].get_fitness();
    uint64_t max_index = 0;
    for (uint64_t i = 1; i < pop.size(); i++)
    {
        if (pop[i].get_fitness() > max)
        {
            max = pop[i].get_fitness();
            max_index = i;
        }
    }
    return max_index;
}

// =================================
// End of Chromosome Implementation
// =================================

int main()
{
    /**
     * @mainpage
     * The Genetic algorithm tries to scan the solution space as much as possible to
     * find an estimate of the optimal solution.
     */

    timer t; // an object helping to calculate the run time
    progressbar<uint64_t> prg;

    // reading the data
    ifstream input("TData.csv");
    if (!input.is_open())
    {
        cout << "\nError opening file!\n";
        return -1;
    }

    // finding number of the columns and saving the headers
    string headers;
    getline(input, headers);
    string s;
    istringstream string_stream(headers);
    uint64_t NCols = 0;                  // number of the columns
    uint64_t NChar = headers.size() + 2; // number of the characters in line 1 (headers) + 2 characters of "\n"
    while (getline(string_stream, s, ','))
    {
        NCols++; // updating the number of columns
    }

    // finding the number of rows
    // the first row was the headers (column names), so we make sure to read the data from the second row
    input.seekg(NChar, ios::beg);
    uint64_t NRows = 0;
    while (getline(input, s))
    {
        NRows++;
    }
    if (input.eof())
        cout << "\nFound the number of rows and columns:\n \t Rows: "
             << NRows << "\tColumns: " << NCols << "\n\n";
    else if (input.fail())
    {
        cout << "Encountered an error in input!\n";
        input.close();
        return -1;
    }
    input.close();
    // updating and printing the progress
    try
    {
        prg.update_progress(10);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }

    // start reading the data again
    ifstream inputt("TData.csv");
    if (!inputt.is_open())
    {
        cout << "\nError opening the file!\n";
        return -1;
    }

    inputt.seekg(NChar, ios::beg); // skipping the headers row

    vector<double> Matrix_elements; // a vector to save the data set
    Matrix_elements.reserve(NRows * (NCols + 1));
    vector<double> nodes(NRows);
    uint64_t i = 0;
    vector<double> v(NCols); // vector to get the output
    s.clear();
    cout << "Started reading the data: ";
    while (getline(inputt, s))
    {
        try
        {
            v = read_D_rows(s, NCols);
            // row number
            nodes[i] = (double)(i + 1);
            Matrix_elements.insert(Matrix_elements.begin() + i * (NCols + 1), (double)(i + 1));

            Matrix_elements.insert(Matrix_elements.begin() + i * (NCols + 1) + 1, v.begin(), v.end());
            // printing the progress (with 10% step size)
            if (i % ((uint64_t)round((double)NRows * 0.1)) == 0)
            {
                cout << "*";
            }
        }
        catch (const exception &e)
        {
            cout << "Error:" << e.what() << '\n';
            inputt.close();
            return -1;
        }
        i++;
    }

    if (inputt.eof())
        cout << "\nReached end of the file.\n\n";
    else if (inputt.fail())
    {
        cout << "Encountered an error in input!\n";
        inputt.close();
        return -1;
    }
    inputt.close();
    // Saving our data set into a matrix
    matrix<double> GeneData(NRows, (NCols + 1), Matrix_elements);
    // updating and printing the progress
    try
    {
        prg.update_progress(25);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }
    prg.print_progress();

    // Solving TSP using the GA
    // calculating the distances between nodes (objective value)
    matrix<double> m_distances = Euclidean_distance(GeneData);
    // updating and printing the progress
    try
    {
        prg.update_progress(30);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }
    prg.print_progress();

    // generating an initial population of size "pop_size"
    uint64_t pop_size = 50;
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
    // updating and printing the progress
    try
    {
        prg.update_progress(35);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }
    prg.print_progress();

    // the algorithm iterates n times
    uint64_t n_iteration = 100;
    // repeating n times
    uint64_t prg5 = 5;
    while (n_iteration--)
    {
        // regenerating another pop_size
        // (50% with cross_over, 20% mutation, 15% simulated anealling, and 15% new random)
        uint64_t pop50 = (uint64_t)round(0.5 * (double)pop_size / 2);
        if (pop50 % 2) // we need pairs of chromosomes to creat children
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
        uint64_t p = 0;       // the index for the population
        uint64_t counter = 0; // a counter of children (new offspring chromosomes)
        try
        {
            while (counter < pop50)
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
            while (counter < pop20)
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
            while (counter < pop15A)
            {
                chromosome<double> child = population[p].annealing(m_distances);
                if (!child.is_in_vector(population))
                {
                    population.push_back(child);
                    p++;
                    counter++;
                }
            }
            counter = 0;
            while (counter < pop15r) // if left
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

            // SELECTION phase
            // calculating the objective value (total distance/cost in TSP) for each individual
            for (chromosome<double> &i : population)
            {
                i.update_fitness(1 / objective_value(i, m_distances));
            }
            // updating and printing the progress every 10 iteration
            if (n_iteration % 10 == 0)
            {
                prg.update_progress(35 + prg5);
                prg.print_progress();
                prg5 += 5;
            }
        }
        // chatching any exception if occurred during the chromosome creation
        catch (const exception &e)
        {
            cout << "\nError: " << e.what() << '\n';
            return -1;
        }

        // choosing 50% of the population with roulette wheel selection mechanism
        double total_objective_values = 0;
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
                    { // has not been added before
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
        population = selected_population;
    }
    // updating and printing the progress
    try
    {
        prg.update_progress(95);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }
    prg.print_progress();

    // now choose the best among the last population (the optimal solution)
    uint64_t max_index = find_max(population);
    // updating and printing the progress
    try
    {
        prg.update_progress(100);
    }
    catch (const exception &e)
    {
        cout << "\nError:" << e.what() << '\n';
        return -1;
    }
    prg.print_progress();
    // printing the estimation of the optimal solution
    cout << "\nOptimal solution is estimated to be:\n"
         << population[max_index] << "with total distance of: "
         << 1 / population[max_index].get_fitness() << " units.\n";
    t.end();
    cout << "\nElapsed time: " << t.seconds() << " seconds.\n\n";
}