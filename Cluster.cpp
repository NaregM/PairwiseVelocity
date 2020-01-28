#include "Cluster.hpp"

vector<double> position(Cluster c)
{
    vector<double> p = {c.x * 1.0e-3, c.y * 1.0e-3, c.z * 1.0e-3};

    return p;
}

vector<double> position_hat(Cluster c)
{
    vector<double> p_hat;
    vector<double> p = {c.x * 1.0e-3, c.y * 1.0e-3, c.z * 1.0e-3};

    double eNorm = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);

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
        s += v1[i] * v2[i];
    }

    return s;
}

istream &operator>>(istream& is, Cluster& coordinates)
{
    is >> coordinates.M >> coordinates.x >> coordinates.y >> coordinates.z >> coordinates.vx >> coordinates.vy >> coordinates.vz;

    return is;
}

void Mcut(vector<Cluster> &C, double Mthr)
{
    auto it = C.begin();
    while (it != C.end())
    {
        if (it->M * 1.0e10 < Mthr)
        {
            it = C.erase(it);

        } else {
            ++it;
        }
    }
}

vector<double> rbins(int nbins, double binsize)
{
    vector<double> rbins_;

    for (int i = 0; i < nbins; i++)
    {
        rbins_.push_back(i*binsize);
    }

    return rbins_;
}

vector<double> dr(Cluster &c1, Cluster &c2)
{
    vector<double> dr_;
    dr_.push_back((c1.x - c2.x) * 1.0e-3);
    dr_.push_back((c1.y - c2.y) * 1.0e-3);
    dr_.push_back((c1.z - c2.z) * 1.0e-3);

    return dr_;
}

double norm(vector<double> v1, vector<double> v2)
{
    return sqrt(dot(v1, v2));
}

vector<double> elementDivie(vector<double> v, double x)
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
