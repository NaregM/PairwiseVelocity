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

vector<Cluster> readFile(string filename, vector<Cluster> Clusters);
double norm(vector<double> &v1, vector<double> &v2);
vector<double> position(Cluster c);
double eDistance(Cluster &c1, Cluster &c2);
vector<double> position_hat(Cluster &c);
double dot(vector<double> v1, vector<double> v2);
istream &operator>>(istream& is, Cluster& coordinates);
void Mcut(vector<Cluster> &C, double Mthr);
void Xadd(vector<Cluster> &C, const double &x);
vector<double> rbins(int nbins, double binsize);
vector<double> dr(Cluster &c1, Cluster &c2);
vector<double> dv(Cluster &c1, Cluster &c2);
vector<double> elementDivid(vector<double> &v, double &x);
vector<double> elementSum(vector<double> &v1, vector<double> &v2);
vector<double> elementSub(vector<double> &v1, vector<double> &v2);
vector<double> elementMultiply(vector<double> &v, double &x);
vector<double> velocity(Cluster &c);
vector<double> v12_direct(vector<Cluster> &clusters, double Mthr, int nbins, double binsize);
vector<double> v12_est_r(vector<Cluster> &clusters, double Mthr, int nbins, double binsize);
vector<double> v12_est_t(vector<Cluster> &clusters, double Mthr, int nbins, double binsize);
void saveResults(string filePath, vector<double> &v);


#endif // CLUSTER
