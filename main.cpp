#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 1000;  //Kuramoto parameters   
double K = 50;                         

double k = 0.1; //Cucker-Smale parameters
double sigma = 1; 
double beta = 1/7;
double R = 10;

double n_c = 10; //Parisi parameter

double t = 0.0;    
size_t nSteps = 1500;
double dt = 0.1;
double T = 5*dt;
double maxdiff = 0.01;

double L = 30;  //box dimension

//**********************//
int Model_type = 1;     ////////  choose model: 0 = Cucker-Smale, metric interaction; 1 = Parisi, topological  ////////  
//**********************//

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
    state_type Int(N);  //interaction terms 

    int n = int(N);

    state_type Phases = Phases_generator(N);
    state_type Xpos = Positions_generator(L, n);
    state_type Ypos = Positions_generator(L, n); 

    std::fstream positions;   //printing positions
    positions.open("Positions.txt", std::ios::out);
    for (int i=0; i<n; i++){   
        positions << Xpos[i] << '\t' << Ypos[i] << '\n';
    }
    positions.close();


    double Adj[n][n];   //Adjacency matrix
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            Adj[i][j] = 0;
        }
    }

    if (Model_type == 0){   //Cucker-Smale
        for (int i=0; i<n; i++){    //filling Adjacency matrix 

            double xi = Xpos[i];
            double yi = Ypos[i];

            for (int j=0; j<n; j++){

                double xj = Xpos[j];
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

                double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                Ri[j] = mod_r;
                if ( j==i ) { Ri[j] = 0; }
            }
            /*std::cout<< "relative distances from " <<i<< " which is in ( " << Xpos[i] << " , " << Ypos[i] << " ) are: " << '\n';
            for (int j=0; j<n; j++){
                std::cout << Ri[j] << '\t';
            }
            std::cout <<'\n';*/

            int i_neighbours = 0;  //conta i vicini di i
            double dr = 0.001;
            double r0 = dr*2;
            double Rmax = L * 1.5;    //"limit" for r; topological interaction has no metric limit, it's only needed in the for loop. L*sqrt(2) would be enough (the box is L*L), 1.5 chosen for certainty

            for (double r=0; r<Rmax; r += r0) {   //NB: dr != r0 needed
                if ( i_neighbours < n_c ){ 
                    
                    for (int j=0; j<n; j++) { 
                        if ( j != i ) {

                            if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                //std::cout << "for " <<i<< " " <<j<< " was seen as a close neighbour as r = " << r <<'\n';
                                i_neighbours += 1;
                                Adj[i][j] = 1/n_c;
                                //if ( i_neighbours == n_c ) { std::cout<< "# of neighbours reached for " <<i<< '\n'; }
                            
                            }
                            
                        }
                    }
                }
                
            }

        }
    }

/*
    double sum=0;
    for (int j=0; j < n; j++) {         //normalizing to have sum_j A_ij = 1 (Markov matrix) NB: when the fireflies will be moving, normalize just for t=0
        for (int i=0; i < n; i++) {
            sum += Adj[i][j];
        }
        for (int i=0; i < n; i++) {
            if( Adj[i][j] != 0 ) { Adj[i][j] /= sum; }
            if( isnan( Adj[i][j] ) ) { Adj[i][j] = 0; }
        }
        sum = 0;
    } 
*/
    
    if ( n <= 10 ) { 
        for (int i=0; i<n; i++){        //printing Adjacency matrix
            for (int j=0; j<n; j++){
                std::cout<< Adj[i][j] << '\t';
            }
            std::cout << '\n';
        } 
    }

    for (int i=0; i<N; i++){   //setting initial values for x[i]
        x[i] = Phases[i];
    }

    std::fstream fout;
    fout.open("solutions.txt", std::ios::out); 
    std::fstream sync;
    sync.open("Synchronization.txt", std::ios::out); 

    //create stepper:
    runge_kutta4<state_type> rk4; 

    for ( int ii=0; ii<nSteps; ++ii ){  //Integration loop

        if (t==0) {             //print initial conditions
            fout << t << '\t';
            for (int i=0; i<N; i++) {
                fout << x[i] << '\t';
            }
            fout << '\n';   
        }

        if ( t >= 0 ) {             //printing m-l/N 
            double m = 0;
            double l = 0;

            for (int i=0; i<N; i++){
                if ( normalizer(x[i]) == 0 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore m+=1 " <<'\n';
                    m += 1;
                }
                if ( normalizer(x[i]) == 1 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore l+=1 " <<'\n';
                    l += 1;
                }
            }
            //std::cout << "m value was: " <<m<< " and l value was: " <<l<< " so m-l/N was: " << (m-l)/N << '\n'; 
            
            sync << t <<'\t'<< (m-l)/N << '\n';

        }
        
        t += dt;    //adjourn current time 
        
        rk4.do_step(MCU(), x , t, dt);     //perform one integration step. Solution x is overwritten

        for (int i=0; i<N; ++i){

            if( Int[i] < 0 ) {         //if the state is incoherent i-th element will stay there longer (aka the state will still be that of x_t)
                x[i] = x_t[i];
                //std::cout << i << " was reset in its t-1 state at time " << t << " as Int[i] was " << Int[i] << '\n';
            }

            Int[i] = 0;   //reset Interaction for the new step

        }

        fout << t << '\t';      //print solution at time t    
        for (int i=0; i<N; i++) {
            fout << x[i] << '\t';
        }
        fout << '\n';   

        if ( ( t!=0 ) && ( t == 0.4 || dmod(t , T , 10) == 0 || dmod(t-(4*dt) , T , 10) == 0 ) ) {     //interazione a t = T, T-dt 
            //std::cout << "interaction at t = " <<t << '\n';

            for (int i=0; i<N; ++i) {
                
                x_t[i] = x[i]; //saving states in x_t[i]

                //std::cout << "evaluating interaction term for " <<i<< '\n';

                for (int j=0; j<N; ++j){

                    if (i != j) {                             //needed for Chi(i,i) = +1, not needed for Adj[i][i] = 0 

                        if( Adj[i][j] != 0 ) {

                            Int[i] += /*(1/N)*/ ( Adj[i][j] * Chi(x[i] , x[j], 0.2) ) ;   //saving interaction terms 
                            /*if (t>0 && t<150){
                                
                                    std::cout <<"Adj["<<i<<"]["<<j<<"]"<< " was: " << Adj[i][j] <<'\n';
                                    std::cout <<"Chi(i,j) was " << Chi(x[i] , x[j], maxdiff) << " as i was " <<x[i]<< " and j was " <<x[j] << " .Their normalized values were i: " <<normalizer(x[i]) << " j: " <<normalizer(x[j]) <<'\n';
                                    std::cout <<"the term added to Int["<<i<<"]"<< " was " << Adj[i][j] * Chi(x[i] , x[j], maxdiff ) <<'\n';
                            } */
                        }
                    }   
                }
                //std::cout << "interaction term for " <<i<< " at time " <<t<< " was " << Int[i] << '\n';
            }
        }

    }
    fout.close();  
    sync.close();

}