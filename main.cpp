#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 1000;  //Kuramoto parameters   
double K = 50;                         

double k = 1;   //Cucker-Smale parameters
double sigma = 1; 
double beta = 8.0 ; 
double R = 5.0;

double n_c = 7.0; //Parisi parameter  (# of topological neighbours)

double t = 0.0;    //time related parameters
size_t nSteps = 1000;
double dt = 0.1; 
double T = 5*dt;

double L = 30.0;  //box dimension

//**********************************************************************************************************************************************
//**********************//
int model_type = 0;     ////////  choose model: 0 = global interaction. 1 = Cucker-Smale, metric interaction. 2 = Parisi, topological.   ////////  
//**********************//
//**********************************************************************************************************************************************

typedef std::vector<double> state_type;

/*struct MCU{             //functor to be passed to do_step for the integration
    void operator() (const state_type &x, state_type &dxdt, const double t) {
        
        for(int i=0; i<N; ++i){    
            dxdt[i] = 1/T; 
        }

    }
};*/

int main(){

    state_type x(N);   // Initial condition, vector of N elements (N ODEs)    

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

    if (model_type == 0) {      //global interaction
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {

                if ( i!=j ) {
                    Adj[i][j] = 1;
                }
                else if ( i == j ) {
                    Adj[i][j] = 0;
                }

            }
        }
    }

    if (model_type == 1){   //Cucker-Smale

        double avg_neighbours=0;
        for (int i=0; i<n; i++){    //filling Adjacency matrix 

            double xi = Xpos[i];
            double yi = Ypos[i];

            std::vector<double> Ri(n);

            for (int j=0; j<n; j++){

                double xj = Xpos[j];
                double yj = Ypos[j];
              
                if( j!=i ) {
                    Adj[i][j] = CS_entries(xi, xj, yi, yj, R, k, sigma, beta); 
                } else if ( j==i ) {
                    Adj[i][j] = 0;
                }

                //voglio controllare in media quante vicine vede ognuna una volta fissato R, perché almeno deve vederne n_c per comparare i due modelli, poi non voglio che ne veda una marea sennò troppo easy
                double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                double mod_r = std::sqrt(r_ij);                 
                Ri[j] = mod_r;
                if ( j==i ) { Ri[j] = 0; }

            }

            double i_neighbours = 0; 
            double dr = 0.001;
            double r0 = dr*2;
            double Rmax = R;    

            for (double r=0; r<Rmax; r += r0) {   
                for (int j=0; j<n; j++) { 
                    if ( j != i ) {
                        if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                            i_neighbours += 1;
                        }
                            
                    }
                }                
            }

            avg_neighbours += (i_neighbours / N);
        }
        std::cout << " # of average neighbours in metric interaction was " << avg_neighbours << '\n';

    } 

    else if (model_type == 2){  //Parisi
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

    if ( n <= 20 ) { 
        for (int i=0; i<n; i++){        //printing Adjacency matrix
            for (int j=0; j<n; j++){
                std::cout<< Adj[i][j] << '\t';
            }
            std::cout << '\n';
        } 
    }

    for (int i=0; i<N; i++){   //setting initial values for x[i]
        x[i] = Phases[i];
        Int[i] = 0;
    }

    std::fstream fout;
    fout.open("solutions.txt", std::ios::out); 
    std::fstream sync;
    sync.open("Synchronization.txt", std::ios::out); 
    std::fstream check;
    check.open("Initial_conditions.txt", std::ios::out);

    //create stepper:
    //runge_kutta4<state_type> rk4; 

    int counter=0;

    for ( int ii=0; ii<nSteps; ++ii ){  //Integration loop

        if ( t >= 0 ) {             //printing m-l/N : m = # of fireflies in 0, l = # of fireflies in 1 
            double m = 0;
            double l = 0;

            for (int i=0; i<N; i++){
                if ( normalizer(x[i] ) == 0 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore m+=1 " <<'\n';
                    m += 1;
                }
                if ( normalizer(x[i] ) == 1 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore l+=1 " <<'\n';
                    l += 1;
                }
            }

            if (t==0) {
                check << m << '\t' << l << '\t' << '\n'; 
            }
            
            sync << t <<'\t'<< (m-l)/N << '\n';

        }

        fout << t << '\t';      //print solution at time t    
        for (int i=0; i<N; i++) {
            fout << x[i] << '\t';
        }
        fout << '\n';          
        
        //qui c'era il pezzo (*) che sta in sync.cpp ora (non serve) 

        //****************************interazione a t = T, t = T-dt***************************************//

        if ( ( t!=0 ) && ( t == 0 +(T-dt) || dmod(t , T , 10) == 0 || dmod(t-(T-dt) , T , 10) == 0 ) ) {     
            //std::cout << "interaction at t = " <<t << '\n';

            for (int i=0; i<n; ++i) {
                //std::cout << "evaluating interaction term for " <<i<< '\n';

                for (int j=0; j<n; ++j){

                    Int[i] += (1/N)* Adj[i][j] * ( tanh(x[j]-x[i]) + (1/N)*delta(x[j]-x[i]) )  ;   //saving interaction terms

                    /*        
                    if (t>2 && t<3){
                        if (Adj[i][j] != 0) {         
                            //std::cout <<"Adj["<<i<<"]["<<j<<"]"<< " was: " << Adj[i][j] <<'\n';
                            //std::cout <<"Chi(i,j) was " << Chi(x[j] - x[i]) << " as xj was " <<x[j] << " and xi was " <<x[i] <<'\n';
                            //std::cout <<"the term added to Int["<<i<<"]"<< " was " << (1/N) * Adj[i][j] * Chi(x[j] - x[i]) <<'\n';
                        }
                    }  */ 
                } 
                //std::cout << "interaction term for " <<i<< " at time " <<t<< " was " << Int[i] << '\n';
            }
        } 

        for (int i=0; i<n; i++) {
            if( Int[i] < -0.00000001 ) {  
                counter += 1;
            }      
        }
        /*
            for (int i=0; i<n; i++) {
                if (t>95 && t<97){
                    if (Int[i] < 0 ) {
                        std::cout << "at time t: " <<t<< " negative interaction term for i = " <<i<< '\n';
                        std::cout << " infact its state was: " << x[i] << " ,while its neighbours were " <<'\n';;
                        for (int j=0; j<n; ++j){
                            if (Adj[i][j] != 0) {
                                std::cout << j << " = " << x[j] <<'\n';
                                if ( x[j] > x[i]+ 0.001 || x[j] < x[i] - 0.001 ) {
                                    std::cout << j << " is the one wich can't synchronize: its state is " << x[j] << " and its neighbours are: " <<'\n';
                                    for (int k=0; k<n; k++) {
                                        if (Adj[j][k] != 0) {
                                            std::cout << x[k] << '\n';
                                        }
                                    }
                                    std::cout << " ---------------------------------------------------------- " <<'\n';
                                }
                            }
                        }
                        std::cout <<'\n';
                    }
                }
            } */
        
        if ( counter > 0 ) { std::cout<< "at time " <<t<< "there were " <<counter<< " negative interaction terms on a total of " <<n<< '\n'; }
        if ( counter >= n-1 ) { 
            std::cout << "******ERROR****** : every interaction term was negative; so every firefly stayed in her state and synchronization was impossible to achieve" <<'\n'; 
            /*std::cout << "trying reset" <<'\n';
            state_type x_new = Phases_generator(n);
            for (int i = 0; i<n; ++i){
                x[i] = x_new[i];          
                Int[i] = 0;
            }*/
        }
        counter = 0;  
        
        t += dt;    //adjourn current time   

        for (int i=0; i<n; ++i){       //without integration      
            if ( (Int[i] > -0.0000001) ) {

                x[i] += 0.2;

            }

            Int[i] = 0;
        }

        for (int i=0; i<n; ++i ){
            if ( x[i] > 1.9999 ) { x[i] = 0; }
        }

    }
    check.close();
    fout.close();  
    sync.close();

}