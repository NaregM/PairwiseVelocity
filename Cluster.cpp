#include "Cluster.hpp"

vector<double> position(Cluster c)
{
    vector<double> p = {c.x * 1.0e-3, c.y * 1.0e-3, c.z * 1.0e-3};

    return p;
}

vector<double> position_hat(Cluster &c)
{
    vector<double> p_hat;
    vector<double> p = {c.x * 1.0e-3, c.y * 1.0e-3, c.z * 1.0e-3};

    double eNorm = norm(p, p);

    for (auto p_ : p)
    {
        p_hat.push_back(p_/eNorm);
    }

    return p_hat;
}

double eDistance(Cluster &c1, Cluster &c2)
{
    vector<double> pos1 = position(c1);
    vector<double> pos2 = position(c2);

    return sqrt(pow(pos1[1] - pos2[1], 2) + pow(pos1[2] - pos2[2], 2) + pow(pos1[3] - pos2[3], 2));
}

double dot(vector<double> v1, vector<double> v2)
{
    if (v1.size() != v2.size())
    {
        cout << "Error" << endl;
    }

    int n = v1.size();

    double s = 0.0;

    for (int i = 0; i < n; i++)
    {
        s += (v1[i] * v2[i]);
    }

    return s;
}

void Mcut(vector<Cluster> &C, double Mthr)
{
    auto it = C.begin();
    while (it != C.end())
    {
        if (it->M * (1.0e10)/0.7 < Mthr)
        {
            it = C.erase(it);

        } else {
            ++it;
        }
    }
}

void Xadd(vector<Cluster> &C, const double &x)
{
    auto it = C.begin();
    while (it != C.end())
    {
        it->x = it->x + x;
        ++it;
    }

}


vector<double> rbins(int nbins, double binsize)
{
    vector<double> rbins_(nbins, 0.0);

    for (int i = 0; i < nbins; ++i)
    {
        rbins_[i] = binsize * (double(i)+0.5);
    }

    return rbins_;
}

vector<double> dr(Cluster &c1, Cluster &c2)
{
    vector<double> dr_;
    dr_.push_back((c1.x - c2.x) * 1.0e-3);
    dr_.push_back((c1.y - c2.y) * 1.0e-3);
    dr_.push_back((c1.z - c2.z) * 1.0e-3);

    for (int k = 0; k < 3; ++k)
    {
        if (dr_[k] > 0.5 * 358.0)
        {
            dr_[k] -= 358.0;
        }
        else if (dr_[k] < -0.5 * 358.0)
        {
            dr_[k] += 358.0;
        }
    }

    return dr_;
}

double norm(vector<double> &v1, vector<double> &v2)
{
    return sqrt(dot(v1, v2));
}

vector<double> elementDivid(vector<double> &v, double &x)
{
    vector<double> v_;
    for (auto vv : v)
    {
        v_.push_back(vv/x);
    }

    return v_;
}

vector<double> dv(Cluster &c1, Cluster &c2)
{
    vector<double> dv_;
    dv_.push_back(c1.vx - c2.vx);
    dv_.push_back(c1.vy - c2.vy);
    dv_.push_back(c1.vz - c2.vz);

    return dv_;
}

vector<double> velocity(Cluster &c)
{
    vector<double> v_;
    v_.push_back(c.vx);
    v_.push_back(c.vy);
    v_.push_back(c.vz);

    return v_;
}

vector<double> elementSum(vector<double> &v1, vector<double> &v2)
{
    vector<double> s;
    for (int i = 0; i < v1.size(); i++)
    {
        s.push_back(v1[i] + v2[i]);
    }

    return s;
}

vector<double> v12_direct(vector<Cluster> &clusters, double Mthr, int nbins, double binsize)
{
    // Remove elements with M < Mthr
    //double Mthr;
    //cout << "Enter mass threshold [M_{solar}]: " << endl;
    //cin >> Mthr;
    Mcut(clusters, Mthr);
    Xadd(clusters, 5e6);
    //cout << "Finial size: " << Clusters.size() << endl;

    int nclus = clusters.size();   // Number of clusters

    cout << "Total numner of clusters: " << nclus << endl;

    //int nbins = 20;  // Number of bins
    //double binsize = 4.0;

    //vector<double> rbins_ = rbins(nbins, binsize);

    vector<int> n_of_r(nbins, 0);
    vector<double> v_of_r(nbins, 0.0);

    for (int i = 0; i < clusters.size(); i++)
    {
        if (i % 100 == 0)

            cout << "Calculating cluster : " << i <<  " --- ";

        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                vector<double> dr_ = dr(clusters[i], clusters[j]);

                vector<double> dv_ = dv(clusters[i], clusters[j]);

                double dr_norm = sqrt(dot(dr_, dr_));
                vector<double> dr_hat = elementDivid(dr_, dr_norm);

                double dv_pair = dot(dv_, dr_hat);

                int ibin = (int)floor(dr_norm/binsize);

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
        //if (n_of_r[i] > 0)
        v_of_r[i] = v_of_r[i]/n_of_r[i];
    }

    return v_of_r;
}

vector<double> v12_estimator(vector<Cluster> &clusters, double Mthr, int nbins, double binsize)
{
    // Remove elements with M < Mthr
    //double Mthr;
    //cout << "Enter mass threshold [M_{solar}]: " << endl;
    //cin >> Mthr;
    Mcut(clusters, Mthr);
    Xadd(clusters, 5e6);
    //cout << "Finial size: " << Clusters.size() << endl;

    int nclus = clusters.size();   // Number of clusters

    cout << "Total numner of clusters: " << nclus << endl;

    //int nbins = 20;  // Number of bins
    //double binsize = 4.0;

    //vector<double> rbins_ = rbins(nbins, binsize);

    vector<double> v12_est(nbins, 0.0);
    vector<double> p_squared(nbins, 0.0);

    for (int i = 0; i < clusters.size(); i++)
    {
        if (i % 100 == 0)

            cout << "Claculating cluster : " << i << " --- ";

        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                vector<double> dr_ = dr(clusters[i], clusters[j]);
                vector<double> r1_hat = position_hat(clusters[i]);
                vector<double> r2_hat = position_hat(clusters[j]);

                double dr_norm = sqrt(dot(dr_, dr_));
                vector<double> dr_hat = elementDivid(dr_, dr_norm);

                vector<double> r1_plus_r2 = elementSum(r1_hat, r2_hat);
                double p_AB = 0.5 * dot(dr_hat, r1_plus_r2);

                double S1 = dot(velocity(clusters[i]), r1_hat);
                double S2 = dot(velocity(clusters[j]), r2_hat);

                int ibin = (int)floor(dr_norm/binsize);

                if (ibin < nbins)
                {
                    v12_est[ibin] += ((S1 - S2) * p_AB);
                    p_squared[ibin] += (p_AB * p_AB);
                }

            } else
                  {
                      continue;
                  }
        }
    }

    for (int i = 0; i < v12_est.size(); i++)
    {
        v12_est[i] = v12_est[i]/p_squared[i];
    }

    return v12_est;
}

void saveResults(string filePath, vector<double> &v)
{
    ofstream output_file(filePath);
    ostream_iterator<double> output_iterator(output_file, "\n");
    copy(v.begin(), v.end(), output_iterator);
}


istream &operator>>(istream& is, Cluster& coordinates)
{
    is >> coordinates.M >> coordinates.x >> coordinates.y >> coordinates.z >> coordinates.vx >> coordinates.vy >> coordinates.vz;

    return is;
}



vector<Cluster> readFile(string filename, vector<Cluster> Clusters)
{
    //char filename[] = "gr_cat_250.txt";
    //vector<Cluster> Clusters;

    ifstream ifs(filename);
    if (ifs)
    {
        copy(istream_iterator<Cluster>(ifs), istream_iterator<Cluster>(), back_inserter(Clusters));
    }

    if(ifs.is_open())
    {
        ifs.close();
    }

    return Clusters;
}
