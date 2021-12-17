#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

// =========
// Interface
// =========

template <typename T>
class matrix
{
public:
    // Constructor to create a zero matrix.
    // First argument: number of rows.
    // Second argument: number of columns.
    matrix(const uint64_t &, const uint64_t &);

    // Constructor to create a diagonal matrix from a vector.
    // Argument: a vector containing the elements on the diagonal.
    // Number of rows and columns is inferred automatically.
    matrix(const vector<T> &);

    // Constructor to create a diagonal matrix from an initializer_list.
    // Argument: an initializer_list containing the elements on the diagonal.
    // Number of rows and columns is inferred automatically.
    matrix(const initializer_list<T> &);

    // Constructor to create a matrix from a vector.
    // First argument: number of rows.
    // Second argument: number of columns.
    // Third argument: a vector containing the elements in row-major order.
    matrix(const uint64_t &, const uint64_t &, const vector<T> &);

    // Constructor to create a matrix from an initializer_list.
    // First argument: number of rows.
    // Second argument: number of columns.
    // Third argument: an initializer_list containing the elements in row-major order.
    matrix(const uint64_t &, const uint64_t &, const initializer_list<T> &);

    // Member function to obtain (but not modify) the number of rows in the matrix.
    uint64_t get_rows() const;

    // Member function to obtain (but not modify) the number of columns in the matrix.
    uint64_t get_cols() const;

    // Overloaded operator () to access matrix elements WITHOUT range checking.
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2.
    // First version: allows modification of the element.
    T &operator()(const uint64_t &, const uint64_t &);

    // Overloaded operator () to access matrix elements WITHOUT range checking.
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2.
    // Second version: does not allow modification of the element.
    const T &operator()(const uint64_t &, const uint64_t &) const;

    // Member function to access matrix elements WITH range checking (throws out_of_range via vector::at).
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2.
    // First version: allows modification of the element.
    T &at(const uint64_t &, const uint64_t &);

    // Member function to access matrix elements WITH range checking (throws out_of_range via vector::at).
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2.
    // Second version: does not allow modification of the element.
    const T &at(const uint64_t &, const uint64_t &) const;

    // Exception to be thrown if the number of rows or columns given to the constructor is zero.
    class zero_size : public invalid_argument
    {
    public:
        zero_size() : invalid_argument("Matrix cannot have zero rows or columns!"){};
    };

    // Exception to be thrown if the vector of elements provided to the constructor is of the wrong size.
    class initializer_wrong_size : public invalid_argument
    {
    public:
        initializer_wrong_size() : invalid_argument("Initializer does not have the expected number of elements!"){};
    };

    // Exception to be thrown if two matrices of different sizes are added or subtracted.
    class incompatible_sizes_add : public invalid_argument
    {
    public:
        incompatible_sizes_add() : invalid_argument("Cannot add or subtract two matrices of different dimensions!"){};
    };

    // Exception to be thrown if two matrices of incompatible sizes are multiplied.
    class incompatible_sizes_multiply : public invalid_argument
    {
    public:
        incompatible_sizes_multiply() : invalid_argument("Two matrices can only be multiplied if the number of columns in the first matrix is equal to the number of rows in the second matrix!"){};
    };

private:
    // The number of rows.
    uint64_t rows = 0;

    // The number of columns.
    uint64_t cols = 0;

    // A vector storing the elements of the matrix in flattened (1-dimensional) form.
    vector<T> elements;
};

// Overloaded binary operator << to easily print out a matrix to a stream.
template <typename T>
ostream &operator<<(ostream &, const matrix<T> &);

// Overloaded binary operator + to add two matrices.
template <typename T>
matrix<T> operator+(const matrix<T> &, const matrix<T> &);

// Overloaded binary operator += to add two matrices and assign the result to the first one.
template <typename T>
matrix<T> operator+=(matrix<T> &, const matrix<T> &);

// Overloaded unary operator - to take the negative of a matrix.
template <typename T>
matrix<T> operator-(const matrix<T> &);

// Overloaded binary operator - to subtract two matrices.
template <typename T>
matrix<T> operator-(const matrix<T> &, const matrix<T> &);

// Overloaded binary operator -= to subtract two matrices and assign the result to the first one.
template <typename T>
matrix<T> operator-=(matrix<T> &, const matrix<T> &);

// Overloaded binary operator * to multiply two matrices.
template <typename T>
matrix<T> operator*(const matrix<T> &, const matrix<T> &);

// Overloaded binary operator * to multiply a scalar on the left and a matrix on the right.
template <typename T>
matrix<T> operator*(const T &, const matrix<T> &);

// Overloaded binary operator * to multiply a matrix on the left and a scalar on the right.
template <typename T>
matrix<T> operator*(const matrix<T> &, const T &);

// ==============
// Implementation
// ==============

template <typename T>
matrix<T>::matrix(const uint64_t &_rows, const uint64_t &_cols)
    : rows(_rows), cols(_cols)
{
    if (rows == 0 or cols == 0)
        throw zero_size();
    elements = vector<T>(rows * cols);
}

template <typename T>
matrix<T>::matrix(const vector<T> &_diagonal)
    : rows(_diagonal.size()), cols(_diagonal.size())
{
    if (rows == 0)
        throw zero_size();
    elements = vector<T>(rows * cols);
    for (uint64_t i = 0; i < rows; i++)
        elements[(cols * i) + i] = _diagonal[i];
}

template <typename T>
matrix<T>::matrix(const initializer_list<T> &_diagonal)
    : matrix(vector<T>(_diagonal)) {}

template <typename T>
matrix<T>::matrix(const uint64_t &_rows, const uint64_t &_cols, const vector<T> &_elements)
    : rows(_rows), cols(_cols), elements(_elements)
{
    if (rows == 0 or cols == 0)
        throw zero_size();
    if (_elements.size() != rows * cols)
        throw initializer_wrong_size();
}

template <typename T>
matrix<T>::matrix(const uint64_t &_rows, const uint64_t &_cols, const initializer_list<T> &_elements)
    : matrix(_rows, _cols, vector<T>(_elements)) {}

template <typename T>
uint64_t matrix<T>::get_rows() const
{
    return rows;
}

template <typename T>
uint64_t matrix<T>::get_cols() const
{
    return cols;
}

template <typename T>
T &matrix<T>::operator()(const uint64_t &row, const uint64_t &col)
{
    return elements[(cols * row) + col];
}

template <typename T>
const T &matrix<T>::operator()(const uint64_t &row, const uint64_t &col) const
{
    return elements[(cols * row) + col];
}

template <typename T>
T &matrix<T>::at(const uint64_t &row, const uint64_t &col)
{
    return elements.at((cols * row) + col);
}

template <typename T>
const T &matrix<T>::at(const uint64_t &row, const uint64_t &col) const
{
    return elements.at((cols * row) + col);
}

template <typename T>
ostream &operator<<(ostream &out, const matrix<T> &m)
{
    out << '\n';
    for (uint64_t i = 0; i < m.get_rows(); i++)
    {
        out << "( ";
        for (uint64_t j = 0; j < m.get_cols(); j++)
            out << m(i, j) << '\t';
        out << ")\n";
    }
    return out;
}

template <typename T>
matrix<T> operator+(const matrix<T> &a, const matrix<T> &b)
{
    if ((a.get_rows() != b.get_rows()) or (a.get_cols() != b.get_cols()))
        throw typename matrix<T>::incompatible_sizes_add();
    matrix<T> c(a.get_rows(), a.get_cols());
    for (uint64_t i = 0; i < a.get_rows(); i++)
        for (uint64_t j = 0; j < a.get_cols(); j++)
            c(i, j) = a(i, j) + b(i, j);
    return c;
}

template <typename T>
matrix<T> operator+=(matrix<T> &a, const matrix<T> &b)
{
    a = a + b;
    return a;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m)
{
    matrix<T> c(m.get_rows(), m.get_cols());
    for (uint64_t i = 0; i < m.get_rows(); i++)
        for (uint64_t j = 0; j < m.get_cols(); j++)
            c(i, j) = -m(i, j);
    return c;
}

template <typename T>
matrix<T> operator-(const matrix<T> &a, const matrix<T> &b)
{
    if ((a.get_rows() != b.get_rows()) or (a.get_cols() != b.get_cols()))
        throw typename matrix<T>::incompatible_sizes_add();
    matrix<T> c(a.get_rows(), a.get_cols());
    for (uint64_t i = 0; i < a.get_rows(); i++)
        for (uint64_t j = 0; j < a.get_cols(); j++)
            c(i, j) = a(i, j) - b(i, j);
    return c;
}

template <typename T>
matrix<T> operator-=(matrix<T> &a, const matrix<T> &b)
{
    a = a - b;
    return a;
}

template <typename T>
matrix<T> operator*(const matrix<T> &a, const matrix<T> &b)
{
    if (a.get_cols() != b.get_rows())
        throw typename matrix<T>::incompatible_sizes_multiply();
    matrix<T> c(a.get_rows(), b.get_cols());
    for (uint64_t i = 0; i < a.get_rows(); i++)
        for (uint64_t j = 0; j < b.get_cols(); j++)
            for (uint64_t k = 0; k < a.get_cols(); k++)
                c(i, j) += a(i, k) * b(k, j);
    return c;
}

template <typename T>
matrix<T> operator*(const T &s, const matrix<T> &m)
{
    matrix<T> c(m.get_rows(), m.get_cols());
    for (uint64_t i = 0; i < m.get_rows(); i++)
        for (uint64_t j = 0; j < m.get_cols(); j++)
            c(i, j) = s * m(i, j);
    return c;
}

template <typename T>
matrix<T> operator*(const matrix<T> &m, const T &s)
{
    return s * m;
}