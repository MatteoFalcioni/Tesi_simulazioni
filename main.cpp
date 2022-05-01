#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 1000;  //Kuramoto parameters
double K = 50;

double k = 0.1; //C-S parameters
double sigma = 1; 
double beta = 3;
double R = 10;

double n_c = 5; //Parisi parameter

double t = 0.0;    
size_t nSteps = 1500;
double dt = 0.1;
double T = 5*dt;

double L = 30;  //box dimension

int Model_type = 0;    //choose model: 0 = Cucker-Smale, metric interaction; 1 = Parisi, topological

typedef std::vector<double> state_type;

struct MCU{
    void operator() (const state_type &x, state_type &dxdt, const double t) {
        
        for(int i=0; i<N; ++i){    
            dxdt[i] = 1/T; 
        }

    }
};


int main(){

    state_type x(N);   // Initial condition, vector of N elements (N ODEs)

    state_type x_t(N);  //needed to save the solution before it's overwritten
    state_type Int(N);

    int n = int(N);

    state_type Phases = Phases_generator(N);
    state_type Xpos = Positions_generator(L, n);
    state_type Ypos = Positions_generator(L, n); 

    /*for (int i=0; i<n; i++){    //printing positions
        std::cout << Xpos[i] << '\t';
    }
    std::cout << '\n';*/


    double Adj[n][n];   //Adjacency matrix

    if (Model_type == 0){   //Cucker-Smale
        for (int i=0; i<n; i++){    //filling Adjacency matrix 
            for (int j=0; j<n; j++){

                double xi = Xpos[i];
                double xj = Xpos[j];
                double yi = Ypos[i];
                double yj = Ypos[j];

              
                if( j!=i ) {
                    Adj[i][j] = CS_entries(xi, xj, yi, yj, R, k, sigma, beta); 
                } else if ( j==i ) {
                    Adj[i][j] = 0;
                }
        
            }
        }
    } 
    else if (Model_type == 1){  //Parisi
        for (int i=0; i<n; i++){
            double xi = Xpos[i];
            double yi = Ypos[i];

            std::vector<double> Ri(n);   //vector of relative distances from i

            for (int j=0; j<n; j++){
                double xj = Xpos[j];
                double yj = Ypos[j];
                double mod_r = 0;

                if ( j != i ) {
                    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                    double mod_r = std::sqrt(r_ij);  //modulo di ri-rj
                }

                Ri[j] = mod_r;
            }

            int i_neighbours = 0;  //conta i vicini di i
            double dr = 0.5;
            double R = L + 1e4;    //"limit" for r, much bigger than the box (topological interaction has no metric limit, it's only needed for the loop)  

            for (int j=0; j<n; j++) {
                if ( j != i ){ 
                    
                    for (double r=0; r<R; r+=dr) {
                        if (i_neighbours < n_c){

                            if ( r > Ri[j]-(dr/10) && r < Ri[j]+(dr/10) ) {

                                i_neighbours += 1;
                                Adj[i][j] = 1/n_c;
                            
                            }
                            r+=dr;
                        }
                    }
                }
                if (j == i){ Adj[i][j] = 0; }
                if (Adj[i][j] < 0.00001 ) { Adj[i][j] = 0; }  //otherwise it won't be 0 but a really small number
            }

        }
    }

    if (Model_type == 0) {
        double sum=0;
        for (int j=0; j < n; j++) {         //normalizing to have sum_j A_ij = 1 (Markov matrix) (only for CS ?)
            for (int i=0; i < n; i++) {
                sum += Adj[i][j];
            }
            for (int i=0; i < n; i++) {
                Adj[i][j] /= sum;
                if ( isnan(Adj[i][j]) ) { Adj[i][j] = 0; }  //dividing very small numbers (0 but with double precision) will result in nan
            }
            sum = 0;
        } 
    }

    /*for (int i=0; i<n; i++){        //printing Adjacency matrix
        for (int j=0; j<n; j++){
            std::cout<< Adj[i][j] << '\t';
        }
        std::cout << '\n';
    }*/

    for (int i=0; i<N; i++){   //setting initial values for x[i]
        x[i] = Phases[i];
    }

    //std::fstream fout;
    //fout.open("solutions.txt", std::ios::out); 
    std::fstream sync;
    sync.open("Synchronization.txt", std::ios::out); 

    //create stepper:
    runge_kutta4<state_type> rk4; 

    for ( int ii=0; ii<nSteps; ++ii ){  //Integration loop

        /*if (t==0) {             //print initial conditions
            fout << t << '\t';
            for (int i=0; i<N; i++) {
                fout << x[i] << '\t';
            }
            fout << '\n';   
        }*/

        if ( t >= 0 ) {             //printing m-l/N 
            double m = 0;
            double l = 0;

            for (int i=0; i<N; i++){
                if ( Normalizer(x[i]) == 0 ) {
                    m += 1;
                }
                if ( Normalizer(x[i]) == 1 ) {
                    l += 1;
                }
            }

            sync << t <<'\t'<< (m-l)/N << '\n';

        }
        
        t += dt;    //adjourn current time        
        
        rk4.do_step(MCU(), x , t, dt);     //perform one integration step. Solution x is overwritten

        for (int i=0; i<N; ++i){

            if( Int[i]<0 ) {         //if the state is incoherent i-th element will stay there longer (aka the state will still be that of x_t)
                x[i] = x_t[i];
                //std::cout << "interaction for " << i << " at time " << t << ";it will stay in its state for one more step" <<'\n';
            }

            Int[i] = 0;   //reset Interaction for the new step

        }

        /*fout << t << '\t';      //print solution at time t    
        for (int i=0; i<N; i++) {
            fout << x[i] << '\t';
        }
        fout << '\n';  */ 

        if ( ( t!=0 ) && ( dmod(t , T , 100) == 0 || dmod(t , T-dt , 100) == 0 ) ) {     //interazione a t = T, T-dt

            for (int i=0; i<N; ++i) {
                
                x_t[i] = x[i]; //saving states in x_t[i]

                for (int j=0; j<N; ++j){  
                    if (i != j) {                             //needed for Chi(i,i) = +1, not needed for Adj[i][i] = 0 

                        Int[i] += (K/N) * Adj[i][j] *  Chi(x[i] , x[j]) ;   //saving interaction terms

                    }                                 
                }

            }
        }

        

    }
    //fout.close();  
    sync.close();

}