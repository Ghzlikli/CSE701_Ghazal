#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <cmath>
//#include "matrix.hpp"

using namespace std;

/**
 * @brief Class of csv
 *  to read an input file in csv format.
 * @tparam T type which is usually a double.
 */
template <typename T>
class csv
{
public:
    /**
     * @brief Constructs a new csv object.
     * Gets the metadata, such as number of rows, number of columns, and header names.
     * @param _file_name Name of the csv file to read.
     */
    csv(const string &);

    /**
     * @brief Exception to be thrown if the file cannot be opened.
     */
    class file_notfound : public invalid_argument
    {
    public:
        file_notfound() : invalid_argument("\nCannot open a file with the given name!\n"){};
    };

    /**
     * @brief Exception to be thrown if getting the file encounters an error.
     */
    class input_failed : public invalid_argument
    {
    public:
        input_failed() : invalid_argument("\nEncountered an error in input!\n"){};
    };

    /**
     * @brief Exception to be thrown if the file contains empty lines.
     */
    class empty_line : public invalid_argument
    {
    public:
        empty_line() : invalid_argument("\nPlease remove the empty line and try again!\n"){};
    };

    /**
     * @brief Exception to be thrown if the number of columns is less than expected.
     */
    class less_column : public out_of_range
    {
    public:
        less_column() : out_of_range("\nExpected less columns!\n"){};
    };

    /**
     * @brief Exception to be thrown if the number of columns is more than expected.
     */
    class more_column : public out_of_range
    {
    public:
        more_column() : out_of_range("\nExpected more columns!\n"){};
    };

    /**
     * @brief Exception to be thrownif a number value is not valid.
     */
    class number_invalid : public invalid_argument
    {
    public:
        number_invalid() : invalid_argument("\nThe number is invalid and cannot be converted!\n"){};
    };

    /**
     * @brief Reads the data set line by line,
     * and returns the dataset in matrix format.
     * @return matrix<double> The received dataset.
     */
    matrix<T> read_data();

    /**
     * @brief Gets the number of columns.
     * @return uint64_t NCols.
     */
    uint64_t get_NCols() const;

    /**
     * @brief Gets the number of rows.
     * @return uint64_t NRows.
     */
    uint64_t get_NRows() const;

    /**
     * @brief Gets the header names.
     * @return uint64_t headers.
     */
    string get_header() const;

    /**
     * @brief Gets the row numbers.
     * @return vector<T> containing row numbers.
     */
    vector<T> get_row_numbers();

private:
    /**
     * @brief Checks whether it is a valid number.
     * @param value_st The string of the value.
     * @return True if it is valid.
     * @return False invalid.
     */
    bool is_valid_num(const string &);

    /**
     * @brief Reads a row of the dataset in csv format.
     * (columns are separated by comma).
     * @param line The row of the data set.
     * @return vector<double> The read row.
     */
    vector<T> read_rows(const string &);

    /**
     * @brief The name of the csv file.
     */
    string datafile;
    /**
     * @brief Number of columns.
     */
    uint64_t NCols = 0;
    /**
     * @brief Header (column names, which is the first line).
     */
    string headers;
    /**
     * @brief Number of rows.
     */
    uint64_t NRows = 0;
    /**
     * @brief Length of the header line.
     */
    uint64_t NChar;
    /**
     * @brief A vector to save the row numbers.
     */
    vector<T> Row_numbers;
    /**
     * @brief A vector to save the data set used in converting it into a matrix.
     */
    vector<T> Matrix_elements;
};

// ==============
// Implementation
// ==============

template <typename T>
csv<T>::csv(const string &_file_name) : datafile(_file_name)
{
    // Reading the data
    ifstream input(datafile);
    if (!input.is_open())
    {
        throw typename csv::file_notfound();
    }
    // Finding number of the columns and saving the headers
    getline(input, headers);
    string line;
    istringstream string_stream(headers);
    NChar = headers.size() + 2; // Number of the characters in line 1 (headers) + 2 characters of "\n"
    while (getline(string_stream, line, ','))
    {
        NCols++; // Updating the number of columns
    }

    // Finding the number of rows
    // The first row was the headers (column names), so we make sure to read the data from the second row
    input.seekg(NChar, ios::beg);
    while (getline(input, line))
    {
        // Stopping on the empty lines
        if (line.size() == 2) // The "\n" at the end of each line is 2 characters in length
        {
            cout << "\nNotice: Line number " << NRows + 1 << " is empty. \n";
            throw typename csv::empty_line();
        }
        NRows++;
    }
    if (input.eof())
    {
        cout << "\nData file is successfully received\n"
             << "\nFound the number of rows and columns:\n \t Rows: "
             << NRows << "\tColumns: " << NCols << "\n\n";
    }
    else if (input.fail())
    {
        input.close();
        throw typename csv::input_failed();
    }
    input.close();
}

template <typename T>
bool csv<T>::is_valid_num(const string &value_st)
{
    // Breaking down the string to characters
    vector<char> c_in_str(value_st.c_str(), value_st.c_str() + value_st.size() + 1);
    uint64_t size = c_in_str.size() - 1;

    uint64_t dot = 0;
    for (uint64_t i = 0; i < size; i++)
    {
        if (c_in_str[i] == ' ')
        {
            continue;
        }
        if (c_in_str[i] == '.') // The number is allowed to have only one decimal point
        {
            dot++;
            if (dot > 1)
            {
                return false;
            }
            continue;
        }
        if (i == 0 && c_in_str[i] == '-') // In case of a negative number, skips the minus sign in the begining only
        {
            continue;
        }
        // The other characters should be only digits
        if (isdigit(c_in_str[i]) == 0)
        {
            return false;
        }
    }
    return true;
}

template <typename T>
vector<T> csv<T>::read_rows(const string &line)
{
    vector<T> v;    // The vector for saving the data
    uint64_t j = 0; // Number of columns read for each row
    istringstream string_stream(line);
    string column;

    // Extracting values of the columns from each row
    while (getline(string_stream, column, ','))
    {
        j++;           // Updating the size
        if (j > NCols) // As soon as it gets more than the limit, it stops
        {
            throw typename csv::less_column();
        }
        // Checking whether it is a number
        if (is_valid_num(column))
        {
            // Adding the number to our vector
            try
            {
                v.push_back(stod(column));
            }
            catch (const exception &e)
            {
                throw typename csv::number_invalid();
            }
        }
        else
        {
            throw typename csv::number_invalid();
        }
    }
    // Checking whether the the number of read values is correct (equals the number of columns)
    if (j < NCols)
    {
        throw typename csv::more_column();
    }
    return v;
}

template <typename T>
matrix<T> csv<T>::read_data()
{
    // To start reading the data again
    ifstream input(datafile);
    if (!input.is_open())
    {
        throw typename csv::file_notfound();
    }

    input.seekg(NChar, ios::beg); // Skipping the headers line
    Matrix_elements.reserve(NRows * (NCols + 1));
    Row_numbers = vector<T>(NRows);
    uint64_t i = 0;
    vector<T> v(NCols); // A vector to get the output
    string line;
    cout << "Started reading the data: ";
    while (getline(input, line))
    {
        v = read_rows(line);
        Row_numbers[i] = (T)(i + 1);
        Matrix_elements.insert(Matrix_elements.begin() + i * (NCols + 1), Row_numbers[i]);
        Matrix_elements.insert(Matrix_elements.begin() + i * (NCols + 1) + 1, v.begin(), v.end());
        // Printing the progress (with 10% step size)
        if (i % ((uint64_t)round((double)NRows * 0.1)) == 0)
        {
            cout << "*";
        }
        i++;
    }
    if (input.eof())
    {
        cout << "\nReached end of the file.\n + All the rows are received successfully.\n\n";
    }
    else if (input.fail())
    {
        input.close();
        throw typename csv::input_failed();
    }
    input.close();
    // Saving the data set into a matrix
    matrix<T> DataSet(NRows, (NCols + 1), Matrix_elements);
    return DataSet;
}

template <typename T>
inline uint64_t csv<T>::get_NCols() const{
    return NCols;
}

template <typename T>
inline uint64_t csv<T>::get_NRows() const{
    return NRows;
}

template <typename T>
inline string csv<T>::get_header() const{
    return headers;
}

template <typename T>
inline vector<T> csv<T>::get_row_numbers(){
    return Row_numbers;
}


// ==========================
// End of CSV Implementation
// ==========================