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

    ifstream ifs(filename);
    if (ifs)
    {
        copy(istream_iterator<Cluster>(ifs), istream_iterator<Cluster>(), back_inserter(Clusters));
    }

    if(ifs.is_open())
    {
        ifs.close();
    }

    // Remove elements with M < Mthr
    double Mthr;
    cout << "Enter mass threshold [M_{solar}]: " << endl;
    cin >> Mthr;
    Mcut(Clusters, Mthr);
    //cout << "Finial size: " << Clusters.size() << endl;

    int nclus = Clusters.size();   // Number of clusters

    cout << "Total numner of clusters: " << nclus << endl;

    int nbins = 20;  // Number of bins
    double binsize = 4.0;

    vector<double> rbins_ = rbins(nbins, binsize);

    vector<int> n_of_r(nbins, 0);
    vector<double> v_of_r(nbins, 0.0);

    for (int i = 0; i < Clusters.size(); i++)
    {
        if (i % 100 == 0)

            cout << "Doing cluster : " << i << endl;

        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                vector<double> dr_ = dr(Clusters[i], Clusters[j]);
                //vector<double> r1 = position_hat(Clusters[i]);
                //vector<double> r2 = position_hat(Clusters[j]);

                vector<double> dv_ = dv(Clusters[i], Clusters[j]);

                double drAbs = sqrt(dot(dr_, dr_));
                vector<double> dr_norm = elementDivie(dr_, drAbs);

                double dv_pair = dot(dv_, dr_norm);

                int ibin = (int)floor(drAbs/binsize);

                //cout << ibin << endl;

                if (ibin < nbins)
                {
                    n_of_r[ibin] += 1;
                    v_of_r[ibin] += dv_pair;
                }

            } else
                  {
                      continue;
                  }
        }
    }

    for (int i = 0; i < v_of_r.size(); i++)
    {
        v_of_r[i] = v_of_r[i]/n_of_r[i];
    }

    for (auto v : v_of_r)
    {
        cout << v << ", ";
    }

    cout << endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << duration.count() * 1e-6 << " [s] " << endl;

    ofstream output_file("results.txt");
    ostream_iterator<double> output_iterator(output_file, "\n");
    copy(v_of_r.begin(), v_of_r.end(), output_iterator);


    return 0;
}
