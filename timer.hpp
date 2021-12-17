#include <chrono>
using namespace std;
class timer
{
public:
void start()
{
start_time = chrono::steady_clock::now();
}
void end()
{
elapsed_time = chrono::steady_clock::now() - start_time;
}
double seconds() const
{
return elapsed_time.count();
}
private:
chrono::time_point<chrono::steady_clock> start_time =
chrono::steady_clock::now();
chrono::duration<double> elapsed_time =
chrono::duration<double>::zero();
};