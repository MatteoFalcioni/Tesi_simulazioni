#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include "sync.h"

int main() {
    double x1 = 0;
    double x2 = 1;
    changeState(x1);
    changeState(x2);
    std::cout << x1 << x2 << '\n';
}