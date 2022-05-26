#ifndef SYNC_H
#define SYNC_H

#include <vector>
#include <random>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <iostream>

int dmod(double t, double T, double res);

int normalizer(double theta);

bool sameRange(double Theta_i, double Theta_j);
double Chi(double theta_i, double theta_j, double maxdiff);

std::vector<double> Phases_generator(int N);
std::vector<double> Positions_generator(double L, int N);

double CS_entries(double xi, double xj, double yi, double yj, double R, double k, double sigma, double beta);

void move (std::vector<double>& pos, double L, double res);

#endif 