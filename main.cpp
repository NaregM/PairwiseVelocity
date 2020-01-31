#include <fstream>
#include <iterator>
#include <chrono>
#include "Cluster.hpp"

using namespace std;
using namespace std::chrono;


int main()
{
    auto start = high_resolution_clock::now();


    char filename[] = "gr_cat_250.txt";
    vector<Cluster> Clusters;

    Clusters = readFile(filename, Clusters);
    /*ifstream ifs(filename);
    if (ifs)
    {
        copy(istream_iterator<Cluster>(ifs), istream_iterator<Cluster>(), back_inserter(Clusters));
    }

    if(ifs.is_open())
    {
        ifs.close();
    }*/

    vector<double> V12_direct, V12_est;
    V12_direct = v12_direct(Clusters, 1.0e14, 20, 4.0);
    V12_est = v12_estimator(Clusters, 1.0e14, 20, 4.0);

    vector<double> r = rbins(20, 4.0);

    /*cout << "Direct: " << endl;
    for (auto v : V12_direct)
    {
        cout << v << ", " << endl;
    }*/

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << duration.count() * 1e-6 << " [s] " << endl;

    saveResults("v12_direct_3e13.txt", V12_direct);
    saveResults("v12_est_3313.txt", V12_est);
    saveResults("rbins.txt", r);

    return 0;
}
