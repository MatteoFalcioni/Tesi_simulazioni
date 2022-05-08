#include "sync.h"

#include <vector>
#include <random>
#include <cmath>

int dmod(double t_, double T_, double res) {

    if ( T_-0.01 >t_ ) { return -1; }    //all the 0.01 are needed to account for double imprecision

    else {
        double t1 = t_+0.01;
        double T1 = T_+0.01;
        int t = trunc(t1*res);
        int T = trunc(T1*res);
        return t % T;
    }

}


int normalizer(double Theta) {
    int theta = trunc(Theta + 0.01);       //tronca theta nel range giusto (2.7 --> 2)
    if (theta % 2 == 0) { theta = 0; }  //range pari spento, range dispari acceso
    else if (theta % 2 != 0) { theta = 1; }
    return theta;
}

bool sameRange(double Theta_i, double Theta_j) {  //Chi definita col normalizer, problema delle lucciole che rimangono indietro e non vengono notate

    int theta_i = normalizer(Theta_i);
    int theta_j = normalizer(Theta_j); 

    int Phase_diff = theta_i - theta_j;

    if (Phase_diff == 0) { return true; } //fasi coerenti (0,0) o (1,1)
    if (Phase_diff != 0) { return false; } //fasi incoerenti (1,0) o (0,1)
}

double Chi(double theta_i, double theta_j, double maxdiff) {  //no pure così è sbagliato. devono SIA essere nelle stesso range SIA differire di meno di maxdiff per dirsi sincronizzate

    bool samerange = sameRange(theta_i, theta_j);

    //std::cout << " initially i and j were i = " <<theta_i<< " j = " <<theta_j<< " so sameRange returned " <<samerange<< '\n';
    while (theta_i > 1) { theta_i -= 1; }
    while (theta_j > 1) { theta_j -= 1; } 
    //std::cout << " after -=1 i and j were i = " <<theta_i<< " j = " <<theta_j<< '\n';

    double phase_diff = theta_i - theta_j;
    if ( phase_diff < 0 ) { phase_diff = -phase_diff; }
    //std::cout << " phase difference was " <<phase_diff<< '\n';

    if ( phase_diff < maxdiff + 0.01 && samerange ) {
        //std::cout << "and therefore Chi returned +1" <<'\n';
        return +1;
    } 
    if ( phase_diff > maxdiff + 0.01 || !samerange ) {
        //std::cout << "and therefore Chi returned -1" <<'\n';
        return -1;
    } 
}

/*double Chi(double Theta_i, double Theta_j) {  //Chi definita col normalizer, problema delle lucciole che rimangono indietro e non vengono notate

    int theta_i = normalizer(Theta_i);
    int theta_j = normalizer(Theta_j); 

    int Phase_diff = theta_i - theta_j;

    if (Phase_diff == 0) { return +1; } //fasi coerenti (0,0) o (1,1)
    if (Phase_diff != 0) { return -1; } //fasi incoerenti (1,0) o (0,1)
}*/

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


    
    

    

