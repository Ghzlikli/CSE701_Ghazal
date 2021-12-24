#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <cmath>
#include <numeric>
#include <iomanip>
//#include "matrix.hpp" This is already included in the main cpp file
#include <random>
#include <algorithm>

using namespace std;

/**
 * @brief Class of Chromosomes.
 * @param genes Chromosomes are created by genes.
 * @param fitness Proportional to objective value; inspired by the fact that
 * in nature, the fittest individuals have a higher chance of survival.
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
     * @param _size The size of the new chromosome.
     */
    chromosome(const uint64_t &);

    /**
     * @brief Constructor to create a new chromosome randomly.
     * @param _genes A vector of genes; the size is inferred automatically.
     */
    chromosome(const vector<T> &);

    /**
     * @brief Returns size of the chromosome.
     * @return Size of the chromosome (uint64_t).
     */
    uint64_t size() const;

    /**
     * @brief Overloaded operator [] to access genes elements.
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n The element number
     * @return T&
     */
    T &operator[](const uint64_t &);

    /**
     * @brief Overloaded operator [] to access genes elements.
     * The indices start from 0, similar to a vector;
     * Second version: does not allow modification of the element.
     * @param n The element number.
     * @return T&
     */
    const T &operator[](const uint64_t &) const;

    /**
     * @brief Access genes elements with range checking (throws out_of_range via vector::at).
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n The element number
     * @return T&
     */
    T &at(const uint64_t &);

    /**
     * @brief Access genes elements with range checking (throws out_of_range via vector::at).
     * The indices start from 0, similar to a vector;
     * First version: allows modification of the element.
     * @param n The element number
     * @return T&
     */
    const T &at(const uint64_t &) const;

    /**
     * @brief Returns a read/write iterator that points to the first element in the genes vector.
     * @return vector<T>::iterator
     */
    vector<T>::iterator begin();

    /**
     * @brief Returns a read/write iterator that points to one past the last element in the genes vector.
     * @return vector<T>::iterator
     */
    vector<T>::iterator end();

    /**
     * @brief Gets the genes.
     * @return vector<T> of genes.
     */
    vector<T> get_genes() const;

    /**
     * @brief Updates the fitness value of the chromosome.
     * @param fit_value The given/new fitness value.
     */
    void update_fitness(const double &);

    /**
     * @brief Gets the fitness value of the chromosome.
     * @return The fitness value (double).
     */
    double get_fitness() const;

    /**
     * @brief Copies a chromosome.
     * @return The copied chromosome<T>.
     */
    chromosome<T> ch_copy();

    /**
     * @brief Genereate a specific number of random values.
     * It is needed for creating a new chromosome with genes in a random order.
     * @param n Number of random values needed.
     * @return vector<double> The vector of n radom values.
     */
    vector<double> assign_rand(const uint64_t &);

    /**
     * @brief Exception to be thrown if two chromosomes of different sizes are cross overed.
     */
    class incompatible_sizes : public invalid_argument
    {
    public:
        incompatible_sizes() : invalid_argument("\nCannot continue due to different dimensions!\n"){};
    };

    /**
     * @brief Exception to be thrown if the distance matrix does not
     *  contain the line between the nodes (genes of the chromosome).
     */
    class no_link_found : public invalid_argument
    {
    public:
        no_link_found() : invalid_argument("\nThe distance matrix does not contain the required link!\n"){};
    };

    /**
     * @brief The single point crossover.
     * @param ch The second parent.
     * @return vector<chromosome<T>> Two new children (offspring chromosomes).
     */
    vector<chromosome<T>> cross_over(chromosome<T> &);

    /**
     * @brief The Swap Mutation:
     * swaps two genes of the children after being through a cross-over.
     * @param ch The second parent.
     * @return vector<chromosome<T>> two mutated children (offspring chromosomes).
     */
    vector<chromosome<T>> mutation(chromosome<T> &);

    /**
     * @brief Checks whether the chromosome is in a vector of chromosomes.
     * @param v_ch The vector of chromosomes.
     * @return TRUE if it is in the container; FALSE otherwise.
     */
    bool is_in_vector(const vector<chromosome<T>> &);

    /**
     * @brief Finds a chromosome from the neighborhood.
     * It can be used for Simulated Anealling.
     * @return * chromosome<T> The new chromosome.
     */
    chromosome<T> neighborhood();

    /**
     * @brief Runs a Simulated Annealing optimization on a given chromosome.
     * @param m_distances The matrix of distances, used for calculating the objective value.
     * @param k_max Number of iterations.
     * @param Temp0 Initial temperature.
     * @param alpha Cool-down rate.
     * @return chromosome<T> The final found chromosome.
     */
    chromosome<T> annealing(matrix<double> const &, const uint64_t &, const double &, const double &);

private:
    /**
     * @brief Size of the chromosome
     */
    uint64_t ch_size = 1;
    /**
     * @brief A vector containing the genes (nodes)
     */
    vector<T> genes;
    /**
     * @brief The fitness of a chromosome to contribute to the next generation. 
     * It is proportional to the objective value.
     */
    double fitness = 0;
};

/**
 * @brief Overloaded binary operator << to easily print out a chromosome object to a stream.
 * @tparam T
 * @param out The out stream.
 * @param ch The chromosome to print.
 * @return ostream&
 */
template <typename T>
ostream &operator<<(ostream &, const chromosome<T> &);

/**
 * @brief Generates a population of new random chromosomes.
 * @tparam T
 * @param genes The genes creating the chromosomes.
 * @param T_number The total number of individuals.
 * @return vector<chromosome<T>> The population.
 */
template <typename T>
vector<chromosome<T>> Pop_Generator(vector<T> const &, const uint64_t &);

/**
 * @brief Calculates the total cost of each solution (individual);
 * In TSP, the total cost is the total distances covered.
 * @tparam T
 * @param ch The solution (chromosome).
 * @param cost The matrix containing the cost of different lines.
 * @return T
 */
template <typename T>
T objective_value(chromosome<T> &, matrix<T> const &);

/**
 * @brief Finds the individual with the maximum value of fitness.
 * @tparam T
 * @param pop The population to search.
 * @return uint64_t The index of the individual with the maximum value.
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
    // The start and end points have to be the fixed
    genes[0] = 1;
    genes[ch_size - 1] = 1;
    vector<T> internal_nodes(ch_size - 2); // Other points
    copy(_genes.begin() + 1, _genes.end(), internal_nodes.begin());
    vector<double> rand_g(ch_size - 2);
    // Assign a random value for each gene
    rand_g = chromosome::assign_rand(ch_size - 2);

    uint64_t max_indx = 0; // Index of the maximum value
    double swap = 0;
    T swap2 = 0;

    // In a descending oreder
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
    // Checking for equal sizes
    if (ch_size != parent_size)
    {
        throw typename chromosome<T>::incompatible_sizes();
    }

    // Randomly selecting a point to cross over
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2);
    uint64_t k = uid(mt);

    // Child #1
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
    // The second part of the first parent is going to be randomly placed into the first child
    vector<T> deleted_part(ch_size - k);
    copy(genes.begin() + k, genes.end(), deleted_part.begin());
    vector<double> rand_g;
    // Assign a random value
    rand_g = assign_rand(deleted_part.size() - 1); // Rxcept the last node (ending point is fixed)
    uint64_t max_indx = 0;                         // Index of the maximum value
    double swap = 0;
    T swap2 = 0;

    // In a descending oreder
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

    // Child #2
    // Copying the parents again
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
    // The second part of the second parent is going to be randomly placed into the second child
    copy(ch.begin() + k, ch.end(), deleted_part.begin());
    // Assign a random value
    rand_g = assign_rand(deleted_part.size() - 1); // Except the last node (ending point is fixed)

    // In a descending oreder
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
    // Copying the first parent again
    chromosome<T> parent(ch_size - 1);
    copy(genes.begin(), genes.end(), parent.begin());
    vector<chromosome<T>> children = parent.cross_over(ch);
    // Randomly selecting two points to swap
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2); // Swapping except the first and last one
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
chromosome<T> chromosome<T>::annealing(matrix<double> const &m_distances,
                                       const uint64_t &k_max, const double &Temp0, const double &alpha)
{
    double Temp = Temp0;
    // Copying the current chromosome
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
            ch0 = new_ch.ch_copy(); // Accepting it
        }
        else
        {
            double P = exp(-(f_new - f0) / Temp);
            // Deciding on whether to accept it or not
            random_device rd;
            mt19937 mt(rd());
            uniform_real_distribution<double> urd(0, 1);
            double r = urd(mt);
            if (P >= r)
            {
                ch0 = new_ch.ch_copy(); // Accepting it
            }
        }
        Temp = Temp0 * pow(alpha, k); // Cooling down the temperature
    }
    return ch0;
}

template <typename T>
chromosome<T> chromosome<T>::neighborhood()
{
    // Copying the current chromosome
    chromosome<T> new_ch(ch_size - 1);
    copy(genes.begin(), genes.end(), new_ch.begin());
    // Randomly selecting two points to swap
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<uint64_t> uid(1, ch_size - 2); // Changing except the first and last points
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
        // Checking whether it is already in the population or not
        if (!new_ch.is_in_vector(ch_Generation))
            ch_Generation[i] = new_ch;
        else
            i--; // Repeat to get a new one
    }
    return ch_Generation;
}

template <typename T>
T objective_value(chromosome<T> &ch, matrix<T> const &cost)
{
    uint64_t size_ch = ch.size();
    uint64_t size_m = cost.get_rows();
    T fitness = 0;
    uint64_t count = 0; // Counts the number of lines added
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
                count++; // Found another link
            }
        }
    }
    // Checking if all the required lines are found
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