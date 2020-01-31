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

    vector<double> V12_direct, V12_est_r, V12_est_t;

    double Mthr = 2.5e13;
    int nbins = 30;
    double binsize = 4.0;

    V12_direct = v12_direct(Clusters, Mthr, nbins, binsize);
    V12_est_r = v12_est_r(Clusters, Mthr, nbins, binsize);
    V12_est_t = v12_est_t(Clusters, Mthr, nbins, binsize);

    vector<double> r = rbins(nbins, binsize);

    /*cout << "Direct: " << endl;
    for (auto v : V12_direct)
    {
        cout << v << ", " << endl;
    }*/

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << duration.count() * 1e-6 << " [s] " << endl;

    saveResults("v12_direct_2.5e13.txt", V12_direct);
    saveResults("v12_estr_2.5e13.txt", V12_est_r);
    saveResults("v12_estt_2.5e13.txt", V12_est_t);
    saveResults("rbins.txt", r);

    return 0;
}
