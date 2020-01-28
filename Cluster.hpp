#ifndef CLUSTER
#define CLUSTER

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <cmath>
#include "Cluster.hpp"

using namespace std;

struct Cluster
{
    int id;
    double M;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

double norm(vector<double> v1, vector<double> v2);
vector<double> position(Cluster);
double eDistance(Cluster &c1, Cluster &c2);
vector<double> position_hat(Cluster c);
double dot(vector<double> v1, vector<double> v2);
istream &operator>>(istream& is, Cluster& coordinates);
void Mcut(vector<Cluster> &C, double Mthr);
vector<double> rbins(int nbins, double binsize);
vector<double> dr(Cluster &c1, Cluster &c2);
vector<double> dv(Cluster &c1, Cluster &c2);
vector<double> elementDivie(vector<double> v, double x);


#endif // CLUSTER
