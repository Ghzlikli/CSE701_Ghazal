#include <iostream>
#include <string>

using namespace std;

/**
 * @brief Class of Progressbar
 * to print the progress in discrete steps.
 * @tparam T
 */
template <typename T>
class progressbar
{
public:
    /**
     * @brief Constructs a new progressbar object,
     * with progress equals to 0.
     */
    progressbar();

    /**
     * @brief Constructs a new progressbar object with a given progress.
     * @param _progress The initial progress.
     */
    progressbar(T const &);

    /**
     * @brief Exception to be thrown if the given progress is not between 0 and 100.
     */
    class invalid_progress : public invalid_argument
    {
    public:
        invalid_progress() : invalid_argument("\nProgress percentage should be between 0 and 100!\n"){};
    };

    /**
     * @brief Updates the current progress.
     * @param _progress The new progress.
     */
    void update_progress(T const &);

    /**
     * @brief Prints the progress bar.
     */
    void print_progress();

private:
    /**
     * @brief The current progress.
     */
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
    // Each 2% has a square or space in the progress bar
    uint64_t filled_sq = (uint64_t)(progress / 2);
    uint64_t spaces = 50 - filled_sq;
    cout << "\n[ ";
    while (filled_sq--)
    {
        cout << "*";
    }
    while (spaces--)
    {
        cout << " ";
    }
    cout << " ] progress: " << progress << "%\n";
}

// ==================================
// End of Progressbar Implementation
// ==================================
