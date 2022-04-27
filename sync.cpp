#include "sync.h"

#include <vector>
#include <random>
#include <cmath>

int dmod(double t_, double T_, double res) {

    if (T_>t_) { return -1; }  //non posso fargli ritornare T perchè magari è 0.2 e poi mi da 0

    else {
        int t = int(t_*res);
        int T = int(T_*res);
        return t % T;
    }

}

int Normalizer(double Theta) {
    int theta = int(Theta);       //tronca theta nel range giusto (2.7 --> 2)
    if (theta % 2 == 0) { theta = 0; }  //range pari spento, range dispari acceso
    else { theta = 1; }
    return theta;
}

double Chi(double Theta_i, double Theta_j) {

    int theta_i = Normalizer(Theta_i);
    int theta_j = Normalizer(Theta_j); 

    int Phase_diff = theta_i - theta_j;

    if (Phase_diff == 0) { return +1; } //fasi coerenti (0,0) o (1,1)
    if (Phase_diff != 0) { return -1; } //fasi incoerenti (1,0) o (0,1)
}

std::vector<double> Phases_generator(int N){
    std::vector<double> Phases(N);

    std::random_device rd;  
    std::mt19937 seed(rd()); 
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for (int i=0; i<N; i++){
        double rand = distribution(seed);
        if (rand <= 0.5){ Phases[i] = 0; } 
        else if (rand > 0.5) { Phases[i] = 1; }
    }
    return Phases;
}

std::vector<double> Positions_generator(double L, int N){
    std::vector<double> Positions(N);

    std::random_device rd;  
    std::mt19937 seed(rd()); 

    std::uniform_real_distribution<double> pos_dist(0.0, L);

    for (int i=0; i<N; i++){
        double Xrand = pos_dist(seed);
        Positions[i] = Xrand;
    }

    return Positions;    
}

double CS_entries(double xi, double xj, double yi, double yj, double R, double k, double sigma, double beta){

    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  //modulo quadro di ri-rj
    double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

    if ( mod_r <= R ) {
        return k / std::pow( (sigma*sigma + r_ij) , beta ) ;
    } else if ( mod_r > R ) {
        return 0;
    }

}


    
    

    
