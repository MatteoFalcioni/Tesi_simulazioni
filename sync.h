#ifndef SYNC_H
#define SYNC_H

#include <vector>
#include <random>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <iostream>

int dmod(double t_, double T_, double res);

int normalizer(double Theta);

double Chi(double Theta_i, double Theta_j);

std::vector<double> Phases_generator(int N);
std::vector<double> Positions_generator(double L, int N);

double CS_entries(double xi, double xj, double yi, double yj, double R, double k, double sigma, double beta);



#endif 